#!/usr/bin/env python

from subprocess import Popen, PIPE
from collections import Counter
import scipy.sparse as sp
import argparse as argp
import pandas as pd
import numpy as np
import subprocess
import traceback
import pysam
import json
import time
import sys
import re



# Calculate score to identify low-complexity reads using DUST algorithm
# (S>2 should be filtered)
def dust(read):
    triplet_counts = Counter(read[i:i+3] for i in range(len(read) - 2))
    l = len(read) - 2
    S = 0
    for count in triplet_counts.values():
        S += count * (count - 1) / 2 / (l - 1)
    return S


def fast_load_alignments_from_bam(hlaBams, filterLowComplexity):
    readNames_to_aligns = {}
    hlaRefIdMappedSet = set()
    hlaRefID_to_totalMappedReads = {}

    # For all HLA*.bam files in directory
    for bam in hlaBams:
        # bamfile = pysam.AlignmentFile(bam, 'rb')
        # sorted_bam_file = bam.replace('.bam', '.sorted.bam')
        # pysam.sort('-o', sorted_bam_file, '-O', "BAM", bam)
        # bamfile.close()
        # pysam.index(sorted_bam_file)
        # with pysam.AlignmentFile(sorted_bam_file, 'rb') as samfile:
        cmdArgs = ['samtools', 'view', '-c', bam]
        result = subprocess.run(cmdArgs, stdout=subprocess.PIPE, encoding='utf-8')
        num_alignments = int(result.stdout.strip())

        cmdArgs = ['samtools', 'view', bam]
        if sys.version[0] == '2':
            pipe = Popen(cmdArgs, stdout=PIPE)
        else:
            pipe = Popen(cmdArgs, stdout=PIPE, encoding='utf8')

        # loop over lines
        ctr = 0
        for line in pipe.stdout:
            # Get read name from field 0, SAM flags from f1, ref id from f2,
            # position from f3, seq from field 9, and tags from field 11
            line = line.strip().split('\t')
            # num_alignments = samfile.count()
            # for read in samfile.fetch():
            ctr += 1
            if ctr % 10000 == 0:
                print(f"  processing alignment {ctr}/{num_alignments}", end='\r')

            readName = line[0]
            readRefId = line[2]
            readPos = int(line[3])
            readCIGAR = line[5]
            readSeq = line[9]
            readTags = line[11:]

            # read_length = read.query_length

            try:
                editDist = [tag for tag in readTags if tag.startswith('NM')][0].split(':')[-1]
            except:
                print('Error parsing tags:')
                print(readName)
                print(readTags)
                print(line)
                print()
                print(traceback.format_exc())
                continue

            # Add this read to the dictionary
            cigarList = list(filter(None, re.split('(\D+)', readCIGAR)))
            alignedSeq = ''
            pos = 0

            for cigar in zip(cigarList[0::2], cigarList[1::2]):
                clen = int(cigar[0])
                if cigar[1] in 'M=XIP':
                    alignedSeq += readSeq[pos:pos + clen]
                    pos += clen
                elif cigar[1] in 'SH':
                    pos += clen

            # if len(alignedSeq) != read_length:
            #     print()

            # Get proper length of matching using corrected readlength
            error_length = int(editDist)
            match_length = len(alignedSeq) - error_length

            if filterLowComplexity:
                passDust = dust(alignedSeq) <= 2
            else:
                passDust = True

            # Disallow clipping on both ends
            if cigarList[1] in 'HS' and cigarList[-1] in 'HS':
                passDust = False

            hlaRefIdMappedSet.add(readRefId)

            if readName not in readNames_to_aligns:
                readNames_to_aligns[readName] = {}
            if readRefId not in readNames_to_aligns[readName]:
                readNames_to_aligns[readName][readRefId] = {
                    'passDust': passDust,
                    'match_length': match_length,
                    'error_length': error_length,
                    'readPos': readPos,
                    'readCIGAR': readCIGAR
                }
            else:
                if match_length > readNames_to_aligns[readName][readRefId]['match_length']:
                    readNames_to_aligns[readName][readRefId] = {
                        'passDust': passDust or readNames_to_aligns[readName][readRefId]['passDust'],
                        'match_length': match_length,
                        'error_length': error_length,
                        'readPos': readPos,
                        'readCIGAR': readCIGAR
                    }

            if readRefId in hlaRefID_to_totalMappedReads:
                hlaRefID_to_totalMappedReads[readRefId] += 1
            else:
                hlaRefID_to_totalMappedReads[readRefId] = 1

    return readNames_to_aligns, hlaRefIdMappedSet, hlaRefID_to_totalMappedReads


def create_hla_read_matrix(readNames_to_aligns):
    read_names = list(sorted(readNames_to_aligns.keys()))
    ref_ids = set()

    # Get list of all ref IDs observed
    for read_align in readNames_to_aligns.values():
        ref_ids.update(read_align.hlaRefID_to_AlignInfo.keys())

    ref_ids = sorted(ref_ids)
    num_reads = len(read_names)
    num_ref_ids = len(ref_ids)

    # Create a mapping from ref_id to index to speed up indexing
    ref_id_to_index = {ref_id: i for i, ref_id in enumerate(ref_ids)}
    read_name_to_index = {read_name: i for i, read_name in enumerate(read_names)}

    match_length_matrix = np.zeros((num_ref_ids, num_reads), dtype=int)
    error_length_matrix = np.zeros((num_ref_ids, num_reads), dtype=int)
    ambig_array = np.zeros(num_reads, dtype=int)
    pass_dust_array = np.zeros(num_reads, dtype=int)

    # matrices have ref_ids as rows and read_names as columns
    # ambiguous and pass_dust arrays have length equal to the number of read_names
    for j, read_name in enumerate(read_names):
        read_align = readNames_to_aligns[read_name]
        read_alignments = read_align.hlaRefID_to_AlignInfo

        ambiguous = read_align.ambig
        pass_dust = any(x for x in read_align.passDust)

        ambig_array[j] = ambiguous
        pass_dust_array[j] = pass_dust

        for ref_id in read_alignments.keys():
            ref_id_index = ref_id_to_index[ref_id]

            alignment = read_alignments[ref_id]
            lm = sum(align.match_length for align in alignment if align)
            em = sum(align.error_length for align in alignment if align)

            match_length_matrix[ref_id_index, j] = lm
            error_length_matrix[ref_id_index, j] = em

    return match_length_matrix, error_length_matrix, ambig_array, pass_dust_array, ref_ids, read_names, ref_id_to_index, read_name_to_index




# Convert read alignment information into matrixes and arrays
def fast_create_hla_read_matrix(readNames_to_aligns):
    read_names = list(sorted(readNames_to_aligns.keys()))
    ref_ids = set()

    # Get list of all ref IDs observed
    for readName, values in readNames_to_aligns.items():
        ref_ids.update(values.keys())

    ref_ids = sorted(ref_ids)
    num_reads = len(read_names)
    num_ref_ids = len(ref_ids)

    # Create a mapping from ref_id to index to speed up indexing
    read_name_to_index = {read_name: i for i, read_name in enumerate(read_names)}
    ref_id_to_index = {ref_id: i for i, ref_id in enumerate(ref_ids)}

    match_length_matrix = np.zeros((num_ref_ids, num_reads), dtype=int)
    error_length_matrix = np.zeros((num_ref_ids, num_reads), dtype=int)
    ambig_array = np.zeros(num_reads, dtype=int)
    pass_dust_array = np.zeros(num_reads, dtype=int)

    # matrices have ref_ids as rows and read_names as columns
    # ambiguous and pass_dust arrays have length equal to the number of read_names
    for j, read_name in enumerate(read_names):
        read_alignments = readNames_to_aligns[read_name]

        for ref_id, alignment_dict in read_alignments.items():
            ref_id_index = ref_id_to_index[ref_id]
            match_length = alignment_dict['match_length']
            error_length = alignment_dict['error_length']
            pass_dust = alignment_dict['passDust']
            pass_dust_array[j] = pass_dust

            match_length_matrix[ref_id_index, j] = match_length
            error_length_matrix[ref_id_index, j] = error_length


    # ambig_array[j] = len(read_alignments['readRefId']) > 1
    # pass_dust_array[j] = any(read_alignments['passDust'])
    ambig_matrix = np.where(match_length_matrix > 0, 1, 0)
    ambig_array = np.sum(ambig_matrix, axis=0) > 1

    match_length_matrix_sparse = sp.csr_matrix(match_length_matrix)
    error_length_matrix_sparse = sp.csr_matrix(error_length_matrix)
    # ambig_array_sparse = sp.csr_matrix(ambig_array.reshape(1, -1))
    # pass_dust_array_sparse = sp.csr_matrix(pass_dust_array.reshape(1, -1))

    # Remove references to dense matrices
    del match_length_matrix
    del error_length_matrix
    # del ambig_array
    # del pass_dust_array

    return match_length_matrix_sparse, error_length_matrix_sparse, ambig_array, pass_dust_array, ref_ids, read_names, ref_id_to_index, read_name_to_index


# Only keep alignments with the maximum match length for each reference

def mask_non_max_values(lm_matrix_sparse, ambig_array):
    # Convert ambig_array_sparse to a column vector

    # Convert lm_matrix_sparse to a dense matrix for max operation
    lm_matrix_dense = lm_matrix_sparse.toarray().astype(float)

    # Find the maximum values along the rows
    max_values = lm_matrix_dense.max(axis=1, keepdims=True)

    # Create a sparse mask for maximum values
    mask = lm_matrix_dense == max_values
    unambig_array = ~ambig_array.astype(bool)
    mask = np.logical_or(mask, unambig_array)

    mask_sparse = sp.csr_matrix(mask)

    del mask
    del lm_matrix_dense
    del max_values

    # Apply the mask to the original sparse matrix
    masked_matrix = sp.csr_matrix(lm_matrix_sparse.multiply(mask_sparse))

    return masked_matrix
# def mask_non_max_values(lm_matrix, ambig_array):
#     match_length_matrix = lm_matrix.copy()
#
#     # max_values = np.max(match_length_matrix, axis=1, keepdims=True)
#     max_values = np.max(match_length_matrix.astype(float), axis=1, keepdims=True)
#
#     mask = match_length_matrix == max_values
#
#     # Also keep all unambiguous alignments
#     unambig_array = ~ambig_array.astype(bool)
#     mask = np.logical_or(mask, unambig_array)
#
#     return np.where(mask, match_length_matrix, 0)


# Replace all non-zero values with the total number of reads mapped to the reference
def createTotalMappedReadsMat(masked_matrix_sparse, ref_ids, hlaRefID_to_totalMappedReads):
    matrix_sparse = masked_matrix_sparse.copy()
    num_rows, num_cols = matrix_sparse.shape

    for i, ref_id in enumerate(ref_ids):
        total_mapped_reads = hlaRefID_to_totalMappedReads[ref_id]

        # Create a dense row vector with the total_mapped_reads value
        row_vector = sp.csr_matrix((np.full(num_cols, total_mapped_reads), (np.zeros(num_cols, dtype=int), np.arange(num_cols))), shape=(1, num_cols))

        # Create a mask for the current row
        row_mask = matrix_sparse[i, :].astype(bool).reshape(1, -1)

        # Apply the total_mapped_reads value where the mask is True
        matrix_sparse[i, :] = row_vector.multiply(row_mask)

    return matrix_sparse
# def createTotalMappedReadsMat(masked_matrix, ref_ids, hlaRefID_to_totalMappedReads):
#     matrix = masked_matrix.copy()
#     for i, ref_id in enumerate(ref_ids):
#         total_mapped_reads = hlaRefID_to_totalMappedReads[ref_id]
#         matrix[i, :] = np.where(matrix[i, :] > 0, total_mapped_reads, 0)
#
#     return matrix


# Keep only alignments with the top n highest number of mapped reads for each reference
def filter_highest_mapped_reads(matrix_sparse, ambig_array, n=1):
    # Convert sparse matrix to a dense matrix
    matrix_dense = matrix_sparse.toarray()

    # Create a boolean mask
    mask = np.zeros_like(matrix_dense, dtype=bool)

    for col_index in range(matrix_dense.shape[1]):
        col = matrix_dense[:, col_index]
        top_n_values = np.unique(np.sort(col)[-n:])
        for value in top_n_values:
            mask[:, col_index] |= (matrix_dense[:, col_index] == value)

    # Also keep all unambiguous alignments
    unambig_array = ~ambig_array.astype(bool)
    mask = np.logical_or(mask, unambig_array)

    # Convert the mask to a sparse matrix
    mask_sparse = sp.csr_matrix(mask)

    del mask
    del matrix_dense

    # Apply the mask to the original sparse matrix
    filtered_matrix_sparse = matrix_sparse.multiply(mask_sparse)

    return filtered_matrix_sparse
# def filter_highest_mapped_reads(matrix, ambig_array, n=1):
#     mat = matrix.copy()
#     mask = np.zeros_like(mat, dtype=bool)
#
#     for col_index in range(mat.shape[1]):
#         col = mat[:, col_index]
#         top_n_values = np.unique(np.sort(col)[-n:])
#         for value in top_n_values:
#             mask[:, col_index] |= (mat[:, col_index] == value)
#
#     # Also keep all unambiguous alignments
#     unambig_array = ~ambig_array.astype(bool)
#     mask = np.logical_or(mask, unambig_array)
#
#     return mat * mask


def load_hla_ref(hlaRefPath):
    hlaRefID_to_seq = {}
    hlaRefID_to_type = {}

    # Read in HLA reference file
    with open(hlaRefPath, 'r') as fHlaRef:
        hlaRef = ''
        refId = ''
        for line in fHlaRef:
            if not line:
                break

            if line[0] == '>':
                if hlaRef:
                    hlaRefID_to_seq[refId] = hlaRef
                    hlaRef = ''
                refId = line.strip().split()[0][1:]
                hlaRefID_to_type[refId] = line.strip().split()[1]

            else:
                hlaRef += line.strip()
        hlaRefID_to_seq[refId] = hlaRef

    return hlaRefID_to_seq, hlaRefID_to_type


def load_alignments_from_bam(hlaBams):
    readNames_to_aligns = {}
    hlaRefIdMappedSet = set()
    hlaRefID_to_totalMappedReads = {}

    # For all HLA*.bam files in directory
    for bam in hlaBams:
        total_fail_dust = 0
        total_pass_dust = 0

        mate = f"{int(bam.split('.')[-4])}"

        # Read the file
        cmdArgs = ['samtools', 'view', bam]
        if sys.version[0] == '2':
            pipe = Popen(cmdArgs, stdout=PIPE)
        else:
            pipe = Popen(cmdArgs, stdout=PIPE, encoding='utf8')

        # loop over lines
        for line in pipe.stdout:
            # Get read name from field 0, SAM flags from f1, ref id from f2,
            # position from f3, seq from field 9, and tags from field 11
            line = line.strip().split('\t')
            [readName, readFlags, readRefId, readPos, readCIGAR, readSeq, readTags] = \
                [line[0], line[1], line[2], int(line[3]), line[5], line[9], line[11:]]
            # read_length = len(readSeq)
            readSeq = readSeq.upper()
            try:
                editDist = [tag for tag in readTags if tag.startswith('NM')][0].split(':')[-1]
            except:
                print('Error parsing tags:')
                print(readName)
                print(readTags)
                print(line)
                print()
                print(traceback.format_exc())
                continue
                # raise
                # sys.exit(1)

            # Add this read to the dictionary
            cigarList = list(filter(None, re.split('(\D+)', readCIGAR)))
            alignedSeq = ''
            pos = 0
            for cigar in zip(cigarList[0::2], cigarList[1::2]):
                clen = int(cigar[0])
                if cigar[1] in 'M=XIP':
                    alignedSeq += readSeq[pos:pos + clen]
                    pos += clen
                elif cigar[1] in 'SH':
                    pos += clen

            # Get proper length of matching using corrected readlength
            error_length = int(editDist)
            match_length = len(alignedSeq) - error_length

            passDust = dust(alignedSeq) <= 2
            # Disallow clipping on both ends
            if cigarList[1] in 'HS' and cigarList[-1] in 'HS':
                passDust = False

            hlaRefIdMappedSet.add(readRefId)
            if readName in readNames_to_aligns:
                readNames_to_aligns[readName].addAlign(readRefId, passDust, mate, match_length, error_length, readPos, readCIGAR)
            else:
                readNames_to_aligns[readName] = ReadAligns(readRefId, passDust, mate, match_length, error_length, readPos, readCIGAR)

            if readRefId in hlaRefID_to_totalMappedReads:
                hlaRefID_to_totalMappedReads[readRefId] += 1
            else:
                hlaRefID_to_totalMappedReads[readRefId] = 1

        while pipe.poll() is None:
            # Process not yet terminated, wait
            time.sleep(0.5)
        if pipe.returncode > 0:
            raise RuntimeError('Error parsing viral-aligned BAM files; aborting.')

    return readNames_to_aligns, hlaRefIdMappedSet, hlaRefID_to_totalMappedReads


def mapReads(hlaBams, hlaRefPath='', annot='', filterLowComplex=False, outputName='hlaType', covMapYmax=0, suppressOutputAndFigures=False):
    # Update coverage arrays
    def update_coverage(refId, readAlign):
        # Get the sequence length for this refId
        seq_len = len(hlaRefID_to_seq[refId])

        # Initialize the coverage array if not already done
        if refId not in hlaRefIdCovArrays:
            hlaRefIdCovArrays[refId] = np.zeros(seq_len, dtype=int)

        # Iterate over the two mates
        for i in range(2):
            mate = readAlign.hlaRefID_to_AlignInfo[refId][i]
            if mate:
                cigarList = list(filter(None, re.split('(\D+)', mate.cigar)))
                pos = mate.pos - 1  # Convert to 0-based indexing

                # Parse the CIGAR string
                opCodes = [op[-1] for op in cigarList[1::2]]
                opLengths = [int(length) for length in cigarList[0::2]]

                # Update the coverage array based on the CIGAR operations
                for opCode, opLength in zip(opCodes, opLengths):
                    if opCode in 'M=X':
                        hlaRefIdCovArrays[refId][pos:pos + opLength] += 1
                        pos += opLength
                    elif opCode in 'DN':
                        pos += opLength

    hlaRefIdGeneDict = {}
    hlaRefIdCovArrays = {}
    hlaRefIdCovDict = {}

    hlaRefID_to_seq, hlaRefID_to_type = load_hla_ref(hlaRefPath)

    readNames_to_aligns, hlaRefIdMappedSet, hlaRefID_to_totalMappedReads = fast_load_alignments_from_bam(hlaBams, filterLowComplex)
    lm_matrix, em_matrix, ambig_array, pass_dust_array, ref_ids, read_names, ref_id_to_index, read_name_to_index = fast_create_hla_read_matrix(readNames_to_aligns)

    if not suppressOutputAndFigures:
        # only create coverage plots for ref_ids mapping to at least a third of the maximum number of reads
        threshold = hlaRefID_to_totalMappedReads.values()
        if len(threshold) > 0:
            threshold = max(threshold) / 3
        else:
            print(hlaRefID_to_totalMappedReads)

        top_ref_ids = {x for x in hlaRefIdMappedSet if hlaRefID_to_totalMappedReads[x] >= threshold}

        for readName, readAlign in readNames_to_aligns.items():
            for refId in readAlign.hlaRefID_to_AlignInfo:
                if refId in top_ref_ids:
                    update_coverage(refId, readAlign)

        for refId, covArray in hlaRefIdCovArrays.items():
            hlaRefIdCovDict[refId] = covArray.tolist()

        with open(f'{outputName}.cov_plot_args.json', 'w') as outFile:
            json_obj = {
                "hlaRefIdCovDict": hlaRefIdCovDict,
                "hlaRefID_to_length": {k: len(v) for k, v in hlaRefID_to_seq.items() if k in top_ref_ids},
                "hlaRefID_to_type": {k: v for k, v in hlaRefID_to_type.items() if k in top_ref_ids},
                "hlaRefIdGeneDict": hlaRefIdGeneDict,
                "outputName": outputName
            }
            json.dump(json_obj, outFile)

    # save lm_matrix and em_matrix as dfs
    # pd.DataFrame(lm_matrix, index=ref_ids, columns=read_names).to_csv(f'{outputName}.lm_matrix.tsv', sep='\t')
    # np.savetxt(f'{outputName}.ambig_array.tsv', ambig_array, fmt='%d')
    # np.savetxt(f'{outputName}.pass_dust_array.tsv', pass_dust_array, fmt='%d')
    # pd.DataFrame(hlaRefID_to_totalMappedReads.items(), columns=['ref_id', 'total_mapped_reads']).to_csv(f"{outputName}.hlaRefID_to_totalMappedReads.tsv", sep='\t')

    # Mask all but the maximum match length for each reference
    masked_mat = mask_non_max_values(lm_matrix, ambig_array)

    # Create a new matrix that replaces all remaining alignments (non-zero values) with the total number of reads mapped to respective reference
    total_reads_mat = createTotalMappedReadsMat(masked_mat, ref_ids, hlaRefID_to_totalMappedReads)

    # Filter out all but the top n read mapping references for each read
    filtered_mat = filter_highest_mapped_reads(total_reads_mat, ambig_array)

    filtered_mat_binary_sparse = sp.csr_matrix((filtered_mat > 0).astype(int))

    if filterLowComplex:
        # Convert pass_dust_array_sparse to a column vector
        pass_dust_array_sparse = pass_dust_array.T

        # Create a mask for low-complexity reads
        fail_dust_mask_dense = ~pass_dust_array_sparse.toarray().astype(bool)
        fail_dust_mask_sparse = sp.csr_matrix(fail_dust_mask_dense.reshape(1, -1))

        # Apply the mask to the binary sparse matrix
        filtered_mat_binary_sparse = filtered_mat_binary_sparse.multiply(~fail_dust_mask_sparse)

    # filtered_mat_binary = np.where(filtered_mat > 0, 1, 0)
    #
    # # Filter out low-complexity reads
    # if filterLowComplex:
    #     fail_dust_mask = ~pass_dust_array.astype(bool)
    #     filtered_mat_binary[:, fail_dust_mask] = 0

    # Create an alignment DataFrame from the filtered matrix

    non_zero_rows = np.flatnonzero(filtered_mat_binary_sparse.sum(axis=1) != 0)
    # non_zero_rows = np.flatnonzero(filtered_mat_binary_sparse.sum(axis=1).toarray() != 0)
    non_zero_cols = np.flatnonzero(filtered_mat_binary_sparse.sum(axis=0) != 0)
    # non_zero_cols = np.flatnonzero(filtered_mat_binary_sparse.sum(axis=0).toarray() != 0)

    # Filter out zero rows and columns from the sparse matrix
    filtered_sparse = filtered_mat_binary_sparse[non_zero_rows, :][:, non_zero_cols]

    # Convert the filtered sparse matrix to a dense NumPy array
    filtered_dense = filtered_sparse.toarray()

    # filter ref_ids and read anmes by non_zero_rows and cols
    ref_ids = [ref_ids[i] for i in non_zero_rows]
    read_names = [read_names[i] for i in non_zero_cols]

    # Create a dense DataFrame from the filtered array
    alignments_dataframe = pd.DataFrame(filtered_dense, index=ref_ids, columns=read_names)

    # df = pd.DataFrame(filtered_mat_binary, index=ref_ids, columns=read_names)
    #
    # # Filter out empty rows and columns
    # alignments_dataframe = df.loc[(df.sum(axis=1) != 0), (df.sum(axis=0) != 0)]

    print()
    print(f"Saving alignment matrix of dimensions {len(alignments_dataframe.index)} references x {len(alignments_dataframe.columns)} reads")
    print()

    mappedCount = 0
    outLine = ''
    nameLine = ''

    # Iterate over the original dictionary
    for readName in alignments_dataframe.columns:
        read_index = read_names.index(readName)

        read_is_ambig = ambig_array[read_index]
        read_passes_dust = pass_dust_array[read_index]

        if not (filterLowComplex and not read_passes_dust):
            mappedCount += 1
            if read_is_ambig:
                outLine += '\t\t\t\tA'
            else:
                outLine += '\t\t\t\tU'
            nameLine += '\t\t\t\t' + readName

    outLine = str(mappedCount) + outLine
    nameLine = nameLine[3:]
    outTable = [nameLine]
    outTable.append(outLine)

    if mappedCount:
        for refId in alignments_dataframe.index:
            hlaName = refId.replace(' ', '')
            outLine = hlaName + " (" + hlaRefID_to_type[refId] + ")"
            for readName in alignments_dataframe.columns:
               if alignments_dataframe.loc[refId, readName] > 0:
                    geneSet = set()
                    ref_id_index = ref_id_to_index[refId]
                    read_name_index = read_name_to_index[readName]

                    match_length = lm_matrix[ref_id_index, read_name_index]
                    error_length = em_matrix[ref_id_index, read_name_index]

                    genes = ','.join(sorted(geneSet))
                    outLine += '\t' + '\t'.join(['1', str(match_length), str(error_length), genes])
               else:
                    outLine += '\t0\t-1\t-1\t'
            outTable.append(outLine)

    with open(outputName + '.mappedReads.tsv', 'w') as outFile:
        for line in outTable:
            outFile.write(str(line) + '\n')

    return outTable


def main(argv):
    # mapParse = argp.ArgumentParser()
    # mapParse.add_argument('bam1')
    # mapParse.add_argument('bam2', nargs='?', help='(optional)', default='not supplied')
    # mapParse.add_argument('-r', '--reference', default=0)
    # mapParse.add_argument('-o', '--outname', type=str, default='./hlaType')
    # mapParse.add_argument('-d', '--disabledust', action='store_true')
    # mapParse.add_argument('-y', '--ylimit', type=int, help='fix a maximum y-value for all coverage map axes', default=0)
    # args = mapParse.parse_args()
    #
    # hlaBams = [args.bam1]
    # if args.bam2 != "not supplied":
    #     hlaBams += [args.bam2]

    hlaBams = ['/Users/zacheliason/Downloads/hla-em/output_paired/trial_0/trial_0.1.Aligned.out.bam']
    ref = '/Users/zacheliason/Downloads/hla-em/hla_gen_ABC.fasta'
    mapReads(hlaBams, hlaRefPath=ref, filterLowComplex=False, outputName='out', suppressOutputAndFigures=True)
    # mapReads(hlaBams, hlaRefPath=args.reference, filterLowComplex=not (args.disabledust), outputName=args.outname, covMapYmax=args.ylimit)


if __name__ == "__main__":
    main(sys.argv)
