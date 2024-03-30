#!/usr/bin/env python

from subprocess import Popen, PIPE
import argparse as argp
import pandas as pd
import numpy as np
import matplotlib
import traceback
import json
import time
import sys
import os
import re


matplotlib.use('Agg')


class AlignInfo:
    def __init__(self, match_length, error_length, pos, cigar):
        self.match_length = match_length
        self.error_length = error_length
        self.pos = pos
        self.cigar = cigar


class ReadAligns:
    def __init__(self, *args):
        if len(args) == 3:
            self.ambig, self.passDust, self.hlaRefID_to_AlignInfo = args

        else:
            refId, passDust, mate, match_length, error_length, pos, cigar = args

            mate = int(mate) - 1
            match_length = int(match_length)
            error_length = int(error_length)
            pos = int(pos)
            self.ambig = False
            self.passDust = [False, False]
            self.passDust[mate] = passDust
            self.hlaRefID_to_AlignInfo = {}
            self.hlaRefID_to_AlignInfo[refId] = [0, 0]
            self.hlaRefID_to_AlignInfo[refId][mate] = AlignInfo(match_length, error_length, pos, cigar)

    def addAlign(self, refId, passDust, mate, match_length, error_length, pos, cigar):
        mate = int(mate) - 1
        match_length = int(match_length)
        error_length = int(error_length)
        pos = int(pos)
        if refId in self.hlaRefID_to_AlignInfo:
            if self.hlaRefID_to_AlignInfo[refId][mate]:
                if match_length > self.hlaRefID_to_AlignInfo[refId][mate].match_length:
                    self.hlaRefID_to_AlignInfo[refId][mate] = AlignInfo(match_length, error_length, pos, cigar)
                    self.passDust[mate] = self.passDust[mate] or passDust
            else:
                self.hlaRefID_to_AlignInfo[refId][mate] = AlignInfo(match_length, error_length, pos, cigar)
                self.passDust[mate] = self.passDust[mate] or passDust
        else:
            self.ambig = True
            self.passDust[mate] = self.passDust[mate] or passDust
            self.hlaRefID_to_AlignInfo[refId] = [0, 0]
            self.hlaRefID_to_AlignInfo[refId][mate] = AlignInfo(match_length, error_length, pos, cigar)


# Calculate score to identify low-complexity reads using DUST algorithm
# (S>2 should be filtered)
def dust(read):
    tripletDict = {}
    for i in range(len(read) - 2):
        c = read[i:i + 3]
        if c in tripletDict:
            tripletDict[c] += 1
        else:
            tripletDict[c] = 1
    S = 0
    l = len(read) - 2
    for trip in tripletDict:
        c = float(tripletDict[trip])
        S += c * (c - 1) / 2 / (l - 1)
    return S


# Convert read alignment information into matrixes and arrays
def create_hla_read_matrix(readNames_to_aligns):
    read_names = list(readNames_to_aligns.keys())
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


# Only keep alignments with the maximum match length for each reference
def mask_non_max_values(lm_matrix, ambig_array):
    match_length_matrix = lm_matrix.copy()

    max_values = np.max(match_length_matrix, axis=1, keepdims=True)
    mask = match_length_matrix == max_values

    # Also keep all unambiguous alignments
    unambig_array = ~ambig_array.astype(bool)
    mask = np.logical_or(mask, unambig_array)

    return np.where(mask, match_length_matrix, 0)


# Replace all non-zero values with the total number of reads mapped to the reference
def createTotalMappedReadsMat(masked_matrix, ref_ids, hlaRefID_to_totalMappedReads):
    matrix = masked_matrix.copy()
    for i, ref_id in enumerate(ref_ids):
        total_mapped_reads = hlaRefID_to_totalMappedReads[ref_id]
        matrix[i, :] = np.where(matrix[i, :] > 0, total_mapped_reads, 0)

    return matrix


# Keep only alignments with the top n highest number of mapped reads for each reference
def filter_highest_mapped_reads(matrix, ambig_array, n=1):
    mat = matrix.copy()
    mask = np.zeros_like(mat, dtype=bool)

    for col_index in range(mat.shape[1]):
        col = mat[:, col_index]
        top_n_values = np.unique(np.sort(col)[-n:])
        for value in top_n_values:
            mask[:, col_index] |= (mat[:, col_index] == value)

    # Also keep all unambiguous alignments
    unambig_array = ~ambig_array.astype(bool)
    mask = np.logical_or(mask, unambig_array)

    return mat * mask


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


def mapReads(hlaBams, hlaRefPath='', annot='', filterLowComplex=True, outputName='hlaType', covMapYmax=0, suppressOutputAndFigures=False):
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

                # # Mark any genes this read covers
                # if refId in hlaRefIdGeneDict:
                #     for gene in hlaRefIdGeneDict[refId]:
                #         gName, gStart, gEnd = gene
                #         gStart, gEnd = int(gStart), int(gEnd)
                #         if np.any(np.logical_and(pos >= gStart, pos < gEnd)):
                #             geneSet.add(gName)

    hlaRefIdGeneDict = {}
    hlaRefIdCovArrays = {}
    hlaRefIdCovDict = {}

    hlaRefID_to_seq, hlaRefID_to_type = load_hla_ref(hlaRefPath)

    readNames_to_aligns, hlaRefIdMappedSet, hlaRefID_to_totalMappedReads = load_alignments_from_bam(hlaBams)

    if not suppressOutputAndFigures:
        # only create coverage plots for ref_ids mapping to at least a fifth of the maximum number of reads
        threshold = max(hlaRefID_to_totalMappedReads.values()) / 3
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

    # Create matrixes and arrays from read alignment information
    lm_matrix, em_matrix, ambig_array, pass_dust_array, ref_ids, read_names, ref_id_to_index, read_name_to_index = create_hla_read_matrix(readNames_to_aligns)

    # Mask all but the maximum match length for each reference
    masked_mat = mask_non_max_values(lm_matrix, ambig_array)

    # Create a new matrix that replaces all remaining alignments (non-zero values) with the total number of reads mapped to respective reference
    total_reads_mat = createTotalMappedReadsMat(masked_mat, ref_ids, hlaRefID_to_totalMappedReads)

    # Filter out all but the top n read mapping references for each read
    filtered_mat = filter_highest_mapped_reads(total_reads_mat, ambig_array)
    filtered_mat_binary = np.where(filtered_mat > 0, 1, 0)

    # Filter out low-complexity reads
    if filterLowComplex:
        fail_dust_mask = ~pass_dust_array.astype(bool)
        filtered_mat_binary[:, fail_dust_mask] = 0

    # Create an alignment DataFrame from the filtered matrix
    df = pd.DataFrame(filtered_mat_binary, index=ref_ids, columns=read_names)

    # Filter out empty rows and columns
    alignments_dataframe = df.loc[(df.sum(axis=1) != 0), (df.sum(axis=0) != 0)]

    print()
    print(f"Saving alignment matrix of dimensions {len(alignments_dataframe.index)} x {len(alignments_dataframe.columns)}")
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


    if not suppressOutputAndFigures:
        with open(outputName + '.mappedReads.tsv', 'w') as outFile:
            for line in outTable:
                outFile.write(str(line) + '\n')

    return outTable


def main(argv):
    mapParse = argp.ArgumentParser()
    mapParse.add_argument('bam1')
    mapParse.add_argument('bam2', nargs='?', help='(optional)', default='not supplied')
    mapParse.add_argument('-r', '--reference', default=0)
    mapParse.add_argument('-o', '--outname', type=str, default='./hlaType')
    mapParse.add_argument('-d', '--disabledust', action='store_true')
    mapParse.add_argument('-y', '--ylimit', type=int, help='fix a maximum y-value for all coverage map axes', default=0)
    args = mapParse.parse_args()

    hlaBams = [args.bam1]
    if args.bam2 != "not supplied":
        hlaBams += [args.bam2]

    # # hlaBams = ['/Users/zacheliason/Downloads/hla-em/output/trial_1/trial_1.1.Aligned.out.bam']
    # hlaBams = ['/Users/zacheliason/Downloads/hla-em/output_paired/trial_0/trial_0.1.Aligned.out.bam']
    # hlaRefPath = '/Users/zacheliason/Downloads/hla_gen.fasta'
    # filterLowComplex = True
    # outname = '/Users/zacheliason/Downloads/hla-em/'
    # outTable = mapReads(hlaBams, hlaRefPath=hlaRefPath, filterLowComplex=filterLowComplex, outputName=outname, suppressOutputAndFigures=False)

    outTable = mapReads(hlaBams, hlaRefPath=args.reference, filterLowComplex=not (args.disabledust), outputName=args.outname, covMapYmax=args.ylimit)


if __name__ == "__main__":
    main(sys.argv)
