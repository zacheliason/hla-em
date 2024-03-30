#!rusr/bin/env python

from src.ManipulateFiles import plot_pie_charts
import pandas as pd
import numpy as np
import matplotlib
import math
import sys
import os
import re


matplotlib.use('Agg')

# Difference value between successive loglikelihoods below which algorithm is deemed to have converged
conVal = 1e-4

numIter = 5
maxSteps = 500


class mappedRead:
    def __init__(self, inputList):
        self.hlaType = inputList[0]
        self.readInfo = []
        self.readNum = 0
        self.genes = []
        count=0
        for val in inputList[1:]:
            if (count%4)==0:
                val = int(val)
                self.readInfo.append([val])
                self.readNum+=val
            elif (count%4)==3:
                if not val:
                    self.genes.append([])
                else:
                    self.genes.append(val.split(','))
            else:
                self.readInfo[-1].append(float(val))
            count+=1

            

def EmAlgo(readsTable, allReadsNum, thresholdTpm=1.5, outputName='hlaType', printResult=True, suppressOutputAndFigures=False):
    iterOut = None
    mappedReads = []
    totalReads = 0
    uniqReads = int(readsTable[1].split('\t')[0])
    isReadAmbig = [x for x in readsTable[1].split('\t')[1:] if x != ""]
    for line in readsTable[2:]:
        line = line.split('\t')
        mappedReads.append(mappedRead(line))
        totalReads += mappedReads[-1].readNum

    if mappedReads:
        # Initialize EM algorithm
        m = len(mappedReads[0].readInfo) #number of reads
        k = len(mappedReads) #number of HLA types

        #Parameters
        for ni in range(numIter):
            if ni==0:
                lOut = -float('inf')
                err = 0.005
                phi = [1./k]*k
            else:
                err = 0.05*np.random.random(1)[0]
                phi = np.random.random(k)
                phi /= phi.sum()

            w=np.zeros([m,k])
            steps=0

            # Calculate initial l
            l = -float('inf')

            converged=False
            while not converged:
                steps+=1
                if(steps>maxSteps):
                    print(f'Iter: {ni}, EM algorithm failed to converge after {maxSteps} steps.')
                    break
                    raise RuntimeError('EM algorithm failed to converge after {} steps; aborting.'.format(maxSteps))
                    
                # E step
                for j,hla in enumerate(mappedReads):
                    for i,readInfo in enumerate(hla.readInfo):
                        if readInfo[0]:
                            Lm = readInfo[1]
                            Le = readInfo[2]
                            w[i,j] = ((1.-err)**Lm * err**Le) * phi[j]
                for i in range(m):
                    w[i,:] = w[i,:]/sum(w[i,:])

                # M step
                ## err
                Bnum = 0
                Bden = 0
                for j,hla in enumerate(mappedReads):
                    for i,readInfo in enumerate(hla.readInfo):
                        if readInfo[0]:
                            Lm = readInfo[1]
                            Le = readInfo[2]
                            Bnum += w[i,j]*Le
                            Bden += w[i,j]*Lm
                B = Bnum/Bden
                err = B/(1.+B)

                ## phi
                for j in range(k):
                    phi[j] = sum(w[:,j])/m

                # Calculate loglikelihood, check change
                l0 = l
                l = 0
                for j,hla in enumerate(mappedReads):
                    for i,readInfo in enumerate(hla.readInfo):
                        if readInfo[0]:
                            Lm = readInfo[1]
                            Le = readInfo[2]
                            try:
                                l +=  w[i,j] * math.log(((1.-err)**Lm * err**Le * phi[j])/w[i,j])
                            except:
                                l +=  w[i,j] * math.log(sys.float_info.min)
                if (l-l0) < conVal:
                    converged=True
                    if ni==0 or l>lOut or iterOut is None:
                        if iterOut is None and ni!=0:
                            print(f'  found {ni} was better')
                        lOut = l
                        errOut = err
                        phiOut = phi
                        stepsOut = steps
                        iterOut = ni

        # Print out results:
        types = []
        typesAll = []
        readProps = []
        emProps = []
        output=[]
        hlaGeneReadCountsDict = {}
        geneNamesSet = set()

        output_lines = []

        # Get number of reads that pass TPM threshold:
        totalOutReads = 0
        for j,hla in enumerate(mappedReads):
            if uniqReads*phiOut[j]*1e6/allReadsNum > thresholdTpm:
                totalOutReads += int(round(uniqReads*phiOut[j]))
            
        for j,hla in enumerate(mappedReads):
            hlaName = hla.hlaType
            if uniqReads*phiOut[j]*1e6/allReadsNum > thresholdTpm:
                types.append(hlaName)
                typesAll.append(hlaName)

                HLA_reference = {
                    'HLAtype': hlaName,
                    'MappedReads': hla.readNum,
                    'MappedProportion': float(hla.readNum) / totalReads,
                    'MLE_Reads': int(round(uniqReads * phiOut[j])),
                    'MLE_Probability': round(uniqReads * phiOut[j])/totalOutReads
                }
                output_lines.append(HLA_reference)

                output.append('{!s}\t{:d}\t{:.5f}\t{:d}\t{:.5f}'.format(hlaName,
                                    hla.readNum, float(hla.readNum)/totalReads,
                                    int(round(uniqReads*phiOut[j])), round(uniqReads*phiOut[j])/totalOutReads))
                readProps.append(float(hla.readNum)/totalReads)

                emProps.append(round(uniqReads*phiOut[j])/totalOutReads)

                # Get per-gene read counts
                if hlaName not in hlaGeneReadCountsDict:
                    hlaGeneReadCountsDict[hlaName] = {}

                for ii in range(len(hla.readInfo)):
                    try:
                        geneList = hla.genes[ii]
                    except:
                        print("Len hla.readInfo: " + str(len(hla.readInfo)))
                        print("current index: " + str(ii))
                        print("error from " + str(hla.genes))
                        sys.exit(1)
                    if isReadAmbig[ii] == 'U':
                        val = 1
                    else:
                        val = phiOut[j]
                    if geneList:
                        for gene in geneList:
                            geneNamesSet.add(gene)
                            if gene not in hlaGeneReadCountsDict[hlaName]:
                                hlaGeneReadCountsDict[hlaName][gene] = val
                            else:
                                hlaGeneReadCountsDict[hlaName][gene] += val
            else:
                if os.path.exists(outputName+'.'+hlaName+'.cov.pdf'):
                    os.remove(outputName+'.'+hlaName+'.cov.pdf')
                if float(hla.readNum)/totalReads > 1e-4:
                    typesAll.append(hlaName)
                    readProps.append(float(hla.readNum)/totalReads)
                    emProps.append(0.0)

        print('Converged to < {:.1e} in {:d} iterations'.format(conVal, stepsOut))

        if output:
            df = pd.DataFrame(output_lines)
            df = df.sort_values('MLE_Probability', ascending=False)
            df.to_csv(outputName + ".results.tsv", sep='\t', index=False)

            # Only print and write results to output file if specified
            if not suppressOutputAndFigures:
                # Write out read counts table
                gene_read_counts_df = pd.DataFrame.from_dict(hlaGeneReadCountsDict, orient='index')
                gene_read_counts_df = gene_read_counts_df.reindex(sorted(gene_read_counts_df.columns), axis=1)
                gene_read_counts_df.fillna(0, inplace=True)
                gene_read_counts_df.to_csv(outputName + '.readCounts.tsv', sep='\t', float_format='%.3f')

            print(f"{steps} steps to converge.")

        else:
            with open(outputName+'.results.tsv','w') as fOut:
                fOut.write('No HLA types detected\n')
            if printResult:
                print('No HLA types detected')
    else:
        with open(outputName+'.results.tsv','w') as fOut:
            fOut.write('No HLA types detected\n')
        if printResult:
            print('No HLA types detected')


def main(argv):
    readTableFile = sys.argv[1]
    readTable = []
    if len(sys.argv) > 3:
        outputName = sys.argv[3]
    else:
        outputName = 'hlaType'
    with open(readTableFile,'r') as inFile:
        for line in inFile:
            readTable.append(line.strip('\n'))
    EmAlgo(readTable, allReadsNum=int(sys.argv[2]), thresholdTpm=1.48, outputName=outputName, printResult=True)

        
if __name__=="__main__":
    main(sys.argv)