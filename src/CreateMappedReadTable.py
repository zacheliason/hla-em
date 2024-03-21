#!/usr/bin/env python
# Create table of mapped read matches and mismatches for use as input to HLA type EM algorithm

from matplotlib import pyplot as plt, lines as lines
from subprocess import Popen, PIPE
import argparse as argp
import pandas as pd
import numpy as np
import matplotlib
import traceback
import time
import sys
import os
import re


matplotlib.use('Agg')

# -o /Users/zacheliason/HLA/Trial18/zachOUT -r /Users/zacheliason/HLA/results/ABC_2_reference.fa /Users/zacheliason/HLA/results/ABC_2MATE.1.Aligned.out.bam /Users/zacheliason/HLA/results/ABC_2MATE.2.Aligned.out.bam
# -o /Users/zacheliason/HLA/Trial18/zachOUT -r /Users/zacheliason/HLA/Trial18/hla_em_fullrun/ABC_2_reference.fa /Users/zacheliason/HLA/Trial18/hla_em_fullrun/ABC_2.1.Aligned.out.bam


class alignInfo:
    def __init__(self, match_length, error_length, pos, cigar):
        self.match_length = match_length
        self.error_length = error_length
        self.pos = pos
        self.cigar = cigar


class readAligns:
    def __init__(self, *args):
        if len(args) == 3:
            self.ambig, self.passDust, self.hlaRefID_to_alignInfo = args

        else:
            refId, passDust, mate, match_length, error_length, pos, cigar = args

            mate = int(mate) - 1
            match_length = int(match_length)
            error_length = int(error_length)
            pos = int(pos)
            self.ambig = False
            self.passDust = [False, False]
            self.passDust[mate] = passDust
            self.hlaRefID_to_alignInfo = {}
            self.hlaRefID_to_alignInfo[refId] = [0, 0]
            self.hlaRefID_to_alignInfo[refId][mate] = alignInfo(match_length, error_length, pos, cigar)

    def addAlign(self, refId, passDust, mate, match_length, error_length, pos, cigar):
        mate = int(mate) - 1
        match_length = int(match_length)
        error_length = int(error_length)
        pos = int(pos)
        if refId in self.hlaRefID_to_alignInfo:
            if self.hlaRefID_to_alignInfo[refId][mate]:
                if match_length > self.hlaRefID_to_alignInfo[refId][mate].match_length:
                    self.hlaRefID_to_alignInfo[refId][mate] = alignInfo(match_length, error_length, pos, cigar)
                    self.passDust[mate] = self.passDust[mate] or passDust
            else:
                self.hlaRefID_to_alignInfo[refId][mate] = alignInfo(match_length, error_length, pos, cigar)
                self.passDust[mate] = self.passDust[mate] or passDust
        else:
            self.ambig = True
            self.passDust[mate] = self.passDust[mate] or passDust
            self.hlaRefID_to_alignInfo[refId] = [0, 0]
            self.hlaRefID_to_alignInfo[refId][mate] = alignInfo(match_length, error_length, pos, cigar)


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


def mapReads(hlaBams, hlaRefPath='', annot='', filterLowComplex=True, outputName='hlaType', covMapYmax=0, suppressOutputAndFigures=True):
    readNames_to_aligns = {}
    hlaRefIdMappedSet = set()
    hlaRefID_to_totalMappedReads = {}
    hlaRefID_to_type = {}
    hlaRefID_to_seq = {}

    hlaRefIdGeneDict = {}
    hlaRefIdCovDict = {}

    installDir = os.path.dirname(os.path.abspath(__file__))
    print("annot: {}".format(annot))
    # Make dict to translate ref seq names (SAM field 2) into HLA type names

    if not annot:
        annot = installDir + '/reference/hla_gene_annot.tsv'

    annotColorDict = {'E1': 'g', 'E2': 'gray', 'E3': 'y', 'E4': 'r', 'E5': 'orange',
                      'E6': 'b', 'E7': 'm', 'E8': 'c', 'L1': 'indigo', 'L2': 'brown'}
    annotColors = ['maroon', 'navy', 'pink', 'g', 'gray', 'k', 'y', 'r', 'orange', 'b', 'm', 'c', 'indigo']

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
    # For all HLA*.bam files in directory

    for bam in hlaBams:
        total_fail_dust = 0
        total_pass_dust = 0
        mate = bam.split('.')[-4]

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
            read_length = len(readSeq)
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
                raise
                sys.exit(1)

            # Add this read to the dictionary
            cigarList = list(filter(None, re.split('(\D+)', readCIGAR)))
            alignedSeq = ''
            pos = 0
            clip_length = 0
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

            if not passDust:
                total_fail_dust += 1
            else:
                total_pass_dust += 1

            hlaRefIdMappedSet.add(readRefId)
            if readName in readNames_to_aligns:
                readNames_to_aligns[readName].addAlign(readRefId, passDust, mate, match_length, error_length, readPos, readCIGAR)
            else:
                readNames_to_aligns[readName] = readAligns(readRefId, passDust, mate, match_length, error_length, readPos, readCIGAR)

            if readRefId in hlaRefID_to_totalMappedReads:
                hlaRefID_to_totalMappedReads[readRefId] += 1
            else:
                hlaRefID_to_totalMappedReads[readRefId] = 1


        while pipe.poll() is None:
            # Process not yet terminated, wait
            time.sleep(0.5)
        if pipe.returncode > 0:
            raise RuntimeError('Error parsing viral-aligned BAM files; aborting.')

    print(total_pass_dust / (total_fail_dust + total_pass_dust))

    # Check if all reads aligned to an HLA reference have equal or better alignment to a reference with more reads
    filteredReadNames_to_aligns = {}
    hlaRefIDs = sorted(hlaRefIdMappedSet)
    for readName in readNames_to_aligns.keys():
        isRedundant = True
        readAlign = readNames_to_aligns[readName]
        for refId in hlaRefIDs:
            if refId in readAlign.hlaRefID_to_alignInfo:
                if not readAlign.ambig:
                    isRedundant = False
                    break

        if isRedundant:
            readAlign_info = readAlign.hlaRefID_to_alignInfo
            try:
                # First we filter out any alignments which do not align to the maximum observed match length

                # initialize array of zeros
                ref_match_lengths = np.zeros(len(readAlign_info))
                ref_Ids = list(map(lambda x: x, readAlign.hlaRefID_to_alignInfo.keys()))

                if len(readAlign_info) > 0:
                    first_key = next(iter(readAlign_info))
                else:
                    continue

                # Add array of match lengths from first (and second, if available) reads
                if readAlign_info[first_key][0]:
                    ref_match_lengths = np.add(ref_match_lengths, np.array(list(map(lambda x: readAlign_info[x][0].match_length if type(readAlign_info[x][0]) == alignInfo else 0, readAlign_info))))
                if readAlign_info[first_key][1]:
                    ref_match_lengths = np.add(ref_match_lengths, np.array(list(map(lambda x: readAlign_info[x][1].match_length if type(readAlign_info[x][1]) == alignInfo else 0, readAlign_info))))

                max_match_length_indices = np.where(ref_match_lengths == np.amax(ref_match_lengths))
                # filter alignments by match length
                refs_max_match_lengths = np.array(ref_Ids)[max_match_length_indices]

                # Now we filter alignments that map to a reference not having the maximum observed reads
                filtered_refs_mapped_nums = {x: hlaRefID_to_totalMappedReads[x] for x in hlaRefID_to_totalMappedReads if x in refs_max_match_lengths}
                ref_mapped_nums = np.array(list(filtered_refs_mapped_nums.values()))
                ref_mapped_nums_ids = np.array(list(filtered_refs_mapped_nums.keys()))

                # Select the n top reference genes that match this read
                n = 2
                # if len(ref_mapped_nums) > n - 1:
                # ind = np.argpartition(ref_mapped_nums, -n)[-n]
                # max_ref_ids = ref_mapped_nums_ids[ref_mapped_nums > ref_mapped_nums[ind] - 1]

                max_mapped_indices = np.where(ref_mapped_nums == np.amax(ref_mapped_nums))
                max_ref_ids = np.array(ref_mapped_nums_ids)[max_mapped_indices]

                hlaRefID_to_alignInfo = {}
                for max_ref_id in max_ref_ids:
                    hlaRefID_to_alignInfo[max_ref_id] = readNames_to_aligns[readName].hlaRefID_to_alignInfo[max_ref_id]

                ambig = readNames_to_aligns[readName].ambig
                passDust = readNames_to_aligns[readName].passDust

                filteredReadNames_to_aligns[readName] = readAligns(ambig, passDust, hlaRefID_to_alignInfo)
            except:
                print(traceback.format_exc())
        else:
            filteredReadNames_to_aligns[readName] = readNames_to_aligns[readName]

    # Now process all reads/read pairs in dict to prepare output table and coverage maps
    # Each column of the outTable is a distinct, mapped read
    # First line of outTable is total mapped read #, followed by unique(U)/ambiguous(A) status of each read/pair
    # There follows 1 line for each HLA reference with at least one read mapped to it.  The first column is the HLA name,
    # and each read column has the following format: [0/1 (whether maps to this reference), match_length (-1 if unmapped), ...
    #  ... error_length (-1 if unmapped), (comma-separated gene list)]
    # First line output
    mappedCount = 0
    outLine = ''
    nameLine = ''
    readNames_to_aligns = filteredReadNames_to_aligns


    filteredReadNames_to_aligns = {}
    filtered_refs = set()

    # Iterate over the original dictionary
    for readName, readAlign in readNames_to_aligns.items():
        # Check conditions for keeping the element
        if not (filterLowComplex and not any(readAlign.passDust)):
            mappedCount += 1
            if readAlign.ambig:
                outLine += '\tA'
            else:
                outLine += '\tU'
            nameLine += '\t' + readName
            # Add the element to the filtered dictionary
            filteredReadNames_to_aligns[readName] = readAlign
            filtered_refs.update(readAlign.hlaRefID_to_alignInfo.keys())

    # Replace the original dictionary with the filtered dictionary
    readNames_to_aligns = filteredReadNames_to_aligns
    outLine = str(mappedCount) + outLine
    outTable = [nameLine]
    outTable.append(outLine)

    # Rest of table
    if mappedCount:
        for refId in filtered_refs:
        # for refId in hlaRefIdMappedSet:
            # print(refId)
            # print(hlaRefID_to_type[refId])
            hlaName = refId.replace(' ', '')
            outLine = hlaName + " (" + hlaRefID_to_type[refId] + ")"
            for readName in readNames_to_aligns:
                readAlign = readNames_to_aligns[readName]
                if refId in readAlign.hlaRefID_to_alignInfo:
                    geneSet = set()
                    match_length = 0
                    error_length = 0

                    # Update read coverage depths for this HLA type
                    if refId not in hlaRefIdCovDict:
                        hlaRefIdCovDict[refId] = [0] * len(hlaRefID_to_seq[refId])
                    for mInd in range(2):
                        mate = readAlign.hlaRefID_to_alignInfo[refId][mInd]
                        if mate:
                            cigarList = list(filter(None, re.split('(\D+)', mate.cigar)))
                            pos = mate.pos
                            for cigar in zip(cigarList[0::2], cigarList[1::2]):
                                if cigar[1] in 'M=X':
                                    for i in range(int(cigar[0])):
                                        # Add to coverage count
                                        try:
                                            hlaRefIdCovDict[refId][pos - 1] += 1
                                        except:
                                            print('readName: {}; refID: {}; startPos: {}; CIGAR: {}; pos: {}'.format(
                                                readName, refId, mate.pos, mate.cigar, pos))
                                            print('error_lengthn(hlaRefIdCovDict[refId]): {}'.format(len(hlaRefIdCovDict[refId])))
                                            raise

                                        # Mark any genes this read covers
                                        if refId in hlaRefIdGeneDict:
                                            for gene in hlaRefIdGeneDict[refId]:
                                                gName = gene[0]
                                                gStart = int(gene[1])
                                                gEnd = int(gene[2])
                                                if gStart <= pos and pos <= gEnd:
                                                    geneSet.add(gName)

                                        pos = pos + 1
                                elif cigar[1] in 'DN':
                                    for i in range(int(cigar[0])):
                                        pos = pos + 1
                                # else :'IPSH'

                            match_length += mate.match_length
                            error_length += mate.error_length
                    genes = ','.join(sorted(geneSet))
                    outLine += '\t' + '\t'.join(['1', str(match_length), str(error_length), genes])
                else:
                    outLine += '\t0\t-1\t-1\t'
            outTable.append(outLine)

    # Plot coverage maps
    if not suppressOutputAndFigures:
        for refId in hlaRefIdCovDict:
            fig = plt.figure(figsize=(9, 4))
            r = fig.canvas.get_renderer()
            cov = fig.add_subplot(111)
            cov.plot(list(range(len(hlaRefIdCovDict[refId]))), hlaRefIdCovDict[refId], 'k', lw=0.8)
            cov.set_ylabel('Read coverage', fontsize=14, color='black')

            hlaName = refId.replace(' ', '')

            plt.title(hlaName + " (" + hlaRefID_to_type[refId] + ")")

            if covMapYmax:
                cov.set_ylim(top=covMapYmax)

            # Plot gene annotations
            glines = []
            glabels = []
            y1end = 0
            y2end = 0
            annotScale = 1.3
            ypos1 = plt.ylim()[0] - (plt.ylim()[1] - plt.ylim()[0]) / 12 * annotScale
            ypos2 = plt.ylim()[0] - (plt.ylim()[1] - plt.ylim()[0]) / 7.9 * annotScale
            ypos3 = plt.ylim()[0] - (plt.ylim()[1] - plt.ylim()[0]) / 5.8 * annotScale
            yposlab1 = plt.ylim()[0] - (plt.ylim()[1] - plt.ylim()[0]) / 8.5 * annotScale
            yposlab2 = plt.ylim()[0] - (plt.ylim()[1] - plt.ylim()[0]) / 6.2 * annotScale
            yposlab3 = plt.ylim()[0] - (plt.ylim()[1] - plt.ylim()[0]) / 4.8 * annotScale
            if refId in hlaRefIdGeneDict:
                ic = 0
                gNameLast = ''
                for gene in hlaRefIdGeneDict[refId]:
                    gName = gene[0]
                    gStart = int(gene[1])
                    gEnd = int(gene[2])

                    tname1 = gName[:2].upper()
                    tname2 = gName[-2:].upper()
                    if (tname1 in annotColorDict and
                            (len(gName) < 3 or gName[2] not in '^*')):
                        gc = annotColorDict[tname1]
                    elif (tname2 in annotColorDict and
                          (len(gName) < 3 or gName[-3] not in '^*')):
                        gc = annotColorDict[tname2]
                    else:
                        if gName != gNameLast:
                            ic += 1
                        gc = annotColors[ic]
                        if ic > 13:
                            ic = 0
                    if gStart >= y1end:
                        ypos = ypos1
                        yposlab = yposlab1
                    elif gStart >= y2end:
                        ypos = ypos2
                        yposlab = yposlab2
                    else:
                        ypos = ypos3
                        yposlab = yposlab3
                    gline = cov.add_line(
                        lines.Line2D([gStart, gEnd], [ypos, ypos], color=gc, clip_on=False, linewidth=2))
                    glines.append(gline)
                    glabel = cov.text(gStart, yposlab, gName)
                    glabels.append(glabel)

                    if ypos == ypos1:
                        y1end = max(gEnd,
                                    cov.transData.inverted().transform(glabel.get_window_extent(renderer=r))[1][0])
                    elif ypos == ypos2:
                        y2end = max(gEnd,
                                    cov.transData.inverted().transform(glabel.get_window_extent(renderer=r))[1][0])
                    gNameLast = gName

            fig.savefig(outputName + '.' + hlaName + '.' + hlaRefID_to_type[refId] + '.cov.pdf', bbox_inches='tight',
                        bbox_extra_artists=glines + glabels)
            plt.close(fig)

    # if not suppressOutputAndFigures:
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

    if (args.reference == 0):
        defaultHlaRef = True
    else:
        defaultHlaRef = False

    outTable = mapReads(hlaBams, hlaRefPath=args.reference, filterLowComplex=not (args.disabledust), outputName=args.outname, covMapYmax=args.ylimit)

    with open(args.outname + '.mappedReads.tsv', 'w') as outFile:
        for line in outTable:
            outFile.write(str(line) + '\n')


if __name__ == "__main__":
    main(sys.argv)
