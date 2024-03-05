#!/usr/bin/env python
import timeit
import argparse as argp
import itertools
import subprocess as subp
import os
import sys
import shutil
import traceback
from whichcraft import which
from EMstep_modified import EmAlgo
from CreateMappedReadTable_modified import mapReads
from TreeInterface import *

__version__ = "1.0.2"


ALL_GENE_GROUPS = "ABC"

def str_to_bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argp.ArgumentTypeError('Boolean value expected.')


def prereqs():
    programs = ["python", "samtools", "STAR"]
    ready = True

    for i in range(0, len(programs)):
        if which(programs[i]) is None:
            print(programs[i] + " not installed. Please install " + programs[i])
            ready = False
    return ready


def cmd(args, write=False, filepath=None):
    if(write == True):
        temp = sys.stdout
        sys.stdout = open(filepath, 'w')

        try:
            subp.check_call(args, stdout=sys.stdout)
        except subp.CalledProcessError as e:
            print("Subprocess error with code: " + str(e.returncode))
            sys.exit(e.returncode)
        except:
            print("An unknown error occurred")
            sys.exit(1)

        sys.stdout = temp

    else:
        try:
            print(' '.join(args))
            subp.check_call(args)
        except subp.CalledProcessError as e:
            print("Subprocess error with code: " + str(e.returncode))
            sys.exit(e.returncode)
        except:
            print("An unknown error occurred")
            exit(1)

    return


def indexReferenceGenes(genomeDir, genomeFastaFiles, genomeSAindexNbases, outname):
    print("Indexing reference genome", flush=True)
    cmd(["STAR",
         "--runMode genomeGenerate",
         "--outTmpDir {}_tmp".format(outname),
         "--genomeDir {path}".format(path=genomeDir),
         "--genomeFastaFiles {path}".format(path=genomeFastaFiles),
         # "--runThreadN {}".format(args.threads),
         "--genomeSAindexNbases {}".format(genomeSAindexNbases)])


def main():
    installDir = os.path.dirname(os.path.abspath(__file__))

    myparse = argp.ArgumentParser(prog='HPV-EM', description='HPV-EM is an HPV genotyping tool that utilizes an expectation maximization algorithm to identify the presence of different HPV genotypes in a sample from RNA-seq data.', formatter_class=lambda prog: argp.RawTextHelpFormatter(prog, width=99999))
    
    # positional arguments
    myparse.add_argument("reads1", help="single-end FASTQ file or first paired-end FASTQ file")
    myparse.add_argument("reads2", nargs='?', help="(optional) second paired-end FASTQ file", default="not supplied")

    # options
    myparse.add_argument('-t','--threads', type=int,  help="number of threads to use [1]", default=1)
    myparse.add_argument('-g','--genomeSAindexNbases', type=int,  help="number of bases to use [6]", default=6)
    myparse.add_argument('-r','--reference', help="viral reference genome in FASTA format,\nto be used in place of default HPV reference", default=0)
    myparse.add_argument('-a','--annotation', help="viral gene annotations in TSV format,\nto be used in place of default HPV annotations\n[{}]".format('$HPV-EMPath/reference/hpv_gene_annot.tsv'), default=installDir+'/reference/hpv_gene_annot.tsv')
    myparse.add_argument('--starviral', help="path to a directory containing STAR-generated\nviral genome indexes based on the above FASTA",default=0)
    myparse.add_argument('--dnaData', type=str_to_bool, help="indicates whether the reads input are RNA or DNA sequences", default=False)
    myparse.add_argument('-o', '--outname', type=str, help="output file name prefix [./hpvEM]", default='./hpvEM')
    myparse.add_argument('-d', '--disabledust', action='store_true', help="disable filtering of low-complexity reads")
    myparse.add_argument('--tpm', type=float, help="TPM threshold for identifying a true positive [1.48]", default=1.48)
    myparse.add_argument('-p', '--printem', action='store_true', help="print EM results to STDOUT")
    myparse.add_argument('-k', '--keepint', action='store_true', help="keep intermediate files")
    myparse.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))

    # other required arguments
    requiredNamed = myparse.add_argument_group('required arguments')
    requiredNamed.add_argument('-s','--stargenome', help="path to a directory containing STAR-generated\nhuman genome indexes", required=True)

    args = myparse.parse_args()

    outPath = args.outname.split('/')
    if len(outPath) > 1:
        outPath = '/'.join(outPath[:-1])
        if not os.path.isdir(outPath):
            cmd(["mkdir", outPath])

    tempGenomeDir = os.path.join(outPath, "temp.fa_STAR")
    if not os.path.isdir(tempGenomeDir):
        os.mkdir(tempGenomeDir)

    tempFasta = os.path.join(outPath, "temp.fa")
    # args.reference, id_to_old = translate_fasta_to_new(args.reference, "/storage1/fs1/jin.zhang/Active/HLA-EM/Ensemble_Genome/zach_HLA_EM/final_dict.json", tempFasta)
    id_to_old = None

    # finding path to reference
    if(args.reference == 0):
        defaultHpvRef = True
        args.reference = installDir + '/reference/combined_pave_hpv.fa'
        args.starviral = installDir + '/reference/combined_pave_hpv_STAR'
    else:
        defaultHpvRef = False
        if args.starviral == 0:
            args.starviral = args.reference + '_STAR'
            if not os.path.isdir(args.starviral):
                indexReferenceGenes(genomeDir=args.starviral, genomeFastaFiles=args.reference, genomeSAindexNbases=args.genomeSAindexNbases, outname=args.outname)
            # print('Please provide the path to a folder of STAR indices based on your specified viral genome using the --starviral argument')
            # sys.exit(1)

    if args.threads < 1:
        args.threads = 1

    if(prereqs() == False):
        sys.exit(1)

    numTiers = 200

    # If the algorithm is using RNA data, skip the last iteration
    if not args.dnaData:
        numTiers = numTiers - 1

    allReadsNum = -1
    hpvBams = []
    argsHumanAlignSTAR = ["STAR",
             "--genomeDir {path}".format(path=args.stargenome),
             "--runThreadN {}".format(args.threads),
             "--chimSegmentMin 18",
             "--outSAMtype BAM Unsorted",
             "--outReadsUnmapped Fastx",
             "--outFilterMultimapNmax 100",
             "--outFilterMismatchNmax 4",
             "--outFileNamePrefix {}.".format(args.outname)]
    if args.reads1.endswith(".gz"):
        argsHumanAlignSTAR.append("--readFilesCommand zcat")

    print("Aligning reads to human genome")
    if args.reads2 == "not supplied":
        cmd(argsHumanAlignSTAR + ["--readFilesIn {}".format(args.reads1)])
    else:
        cmd(argsHumanAlignSTAR + ["--readFilesIn {} {}".format(args.reads1, args.reads2)])

    referenceTree = build_tree_from_fasta(args.reference, tempFasta)
    newick_str = referenceTree.to_newick(referenceTree.get_root())

    try:
        with open(os.path.join(outPath, "newick.nwk"), "w") as n:
            n.write(newick_str)
    except:
        print("NEWICK PRINT FAILED")
        pass

    currentTier = 0
    geneGroup = ALL_GENE_GROUPS
    finalTier = False
    for i in range(currentTier, numTiers):
        print(f"Beginning tier {i + 1} of {numTiers}")

        runString = f"Iteration {geneGroup.lower()}-{i + 1}"
        print(runString)

        if i == numTiers - 1:
            finalTier = True

        suppressOutputAndFigures = True

        print(f"num genes before filtering by breadth: {referenceTree.get_num_genes()}")
        # Filter tree (use i + 1 so that it accurately reflects current search tier)
        num_genes = filter_fasta_by_breadth(referenceTree, tempFasta, i + 1, outpath=outPath, keep_int=True)
        print(f"num genes after filtering by breadth: {num_genes}")

        # TODO find out if indexing should be done only on the first iteration
        # if i > 0:
        # Remove all files in tempGenomeDir as per the instructions found in https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf (Sec. 2)
        for file in os.listdir(tempGenomeDir):
            os.remove(os.path.join(tempGenomeDir, file))
        indexReferenceGenes(genomeDir=tempGenomeDir, genomeFastaFiles=tempFasta,genomeSAindexNbases=args.genomeSAindexNbases, outname=args.outname)
        # Set indexed reference directory to the new temporary directory
        args.starviral = tempGenomeDir

        if args.reads2 == "not supplied":
            with open('{}.Log.final.out'.format(args.outname), 'r') as logFile:
                for line in logFile:
                    line = line.strip()
                    if line.startswith('Number of input reads'):
                        allReadsNum = int(line.split()[-1])
                        print("Total reads: {}".format(allReadsNum))
                        break

                print("Aligning reads to HPV genomes", flush=True)
                cmd(["STAR",
                     "--genomeDir {path}".format(path=args.starviral),
                     "--readFilesIn {sampleName}.Unmapped.out.mate1".format(sampleName=args.outname),
                     "--runThreadN {}".format(args.threads),
                     "--twopassMode Basic",
                     "--outSAMtype BAM Unsorted",
                     "--outSAMattributes NH HI NM MD AS XS",
                     "--outFilterScoreMinOverLread 0",
                     "--outFilterMatchNminOverLread 0",
                     "--outFilterMatchNmin 0",
                     "--outFilterMultimapNmax 999",
                     "--outFilterMismatchNoverLmax 0.08",
                     "--winAnchorMultimapNmax 1000",
                     "--outFileNamePrefix {}.1.".format(args.outname)])
                hpvBams.append('{}.1.Aligned.out.bam'.format(args.outname))

        else:
            # Paired reads block
            with open('{}.Log.final.out'.format(args.outname), 'r') as logFile:
                for line in logFile:
                    line = line.strip()
                    if line.startswith('Number of input reads'):
                        allReadsNum = int(line.split()[-1])
                        break

            print("Aligning reads to HPV genomes", flush=True)
            cmd(["STAR",
                 "--genomeDir {path}".format(path=args.starviral),
                 "--readFilesIn {sampleName}.Unmapped.out.mate1".format(sampleName=args.outname),
                 "--runThreadN {}".format(args.threads),
                 "--twopassMode Basic",
                 "--outSAMtype BAM Unsorted",
                 "--outSAMattributes NH HI NM MD AS XS",
                 "--outFilterMultimapNmax 999",
                 "--outFilterMismatchNmax 999",
                 "--outFilterMismatchNoverLmax 0.08",
                 "--outFileNamePrefix {}.1.".format(args.outname)])
            hpvBams.append('{}.1.Aligned.out.bam'.format(args.outname))
            cmd(["STAR",
                 "--genomeDir {path}".format(path=args.starviral),
                 "--readFilesIn {sampleName}.Unmapped.out.mate2".format(sampleName=args.outname),
                 "--runThreadN {}".format(args.threads),
                 "--twopassMode Basic",
                 "--outSAMtype BAM Unsorted",
                 "--outSAMattributes NH HI NM MD AS XS",
                 "--outFilterMultimapNmax 999",
                 "--outFilterMismatchNmax 999",
                 "--outFilterMismatchNoverLmax 0.08",
                 "--outFileNamePrefix {}.2.".format(args.outname)])
            hpvBams.append('{}.2.Aligned.out.bam'.format(args.outname))


        # TODO comment back in to save bam files
        destBam = os.path.join(outPath, f"{geneGroup}_{str(i + 1)}MATE.1.Aligned.out.bam")
        while os.path.exists(destBam):
            print(f"{destBam} already exists!")
            head, tail = os.path.split(destBam)
            tail = "copy_" + tail
            destBam = os.path.join(head, tail)
        if os.path.exists(os.path.join(outPath, "v1.1.Aligned.out.bam")):
            shutil.copyfile(os.path.join(outPath, "v1.1.Aligned.out.bam"), destBam)
        if os.path.exists(os.path.join(outPath, "v1.2.Aligned.out.bam")):
            shutil.copyfile(os.path.join(outPath, "v1.2.Aligned.out.bam"), destBam.replace("MATE.1", "MATE.2"))
        if os.path.exists(os.path.join(outPath, "temp.fa")):
            shutil.copyfile(os.path.join(outPath, "temp.fa"), os.path.join(outPath, f"{geneGroup}_{str(i + 1)}_reference.fa"))

        # Only output figures and results if on the final iteration
        if finalTier:
            suppressOutputAndFigures = False

        print(f"Creating read table", flush=True)
        readsTable = mapReads(hpvBams, defaultHpvRef=defaultHpvRef, hpvRefPath=tempFasta,
                              filterLowComplex=not (args.disabledust), outputName=args.outname,
                              annot=args.annotation,
                              suppressOutputAndFigures=suppressOutputAndFigures)

        print(f"running EM algorithm", flush=True)
        emResults = EmAlgo(readsTable, allReadsNum, thresholdTpm=args.tpm, outputName=args.outname,
                           printResult=args.printem, suppressOutputAndFigures=suppressOutputAndFigures)


        print(output_results(emResults, i + 1, capture_gene_group=geneGroup, id_to_old=id_to_old))

        # If final iteration, stop iterating
        if finalTier:
            print(f"==={geneGroup}=final=tier===================")
            print("FINAL RESULTS", flush=True)
            print(sort_results_hla(emResults))

            break

        # If not the final iteration, filter the tree using results from the EMA algorithm and continue to the next run
        else:
            referenceTree = build_tree_from_results(emResults, referenceTree, i + 1, tempFasta, final_tier=finalTier)

        # If the number of genes contained in the tree equal the number of emResults, break out of loop
        if not finalTier and referenceTree.get_num_genes() <= referenceTree.get_num_filters():
            print(f"num genes in tree: {referenceTree.get_num_genes()}, num filters used to get them: {referenceTree.get_num_filters()}")
            print(referenceTree.get_num_filters())
            print("===exiting=early===================")
            suppressOutputAndFigures = False
            # TODO (run this logic by Matt and see if it's sound) Rerun the EmAlgorithm, this time with figures and output
            # readsTable = mapReads(hpvBams, defaultHpvRef=defaultHpvRef, hpvRefPath=tempFasta, filterLowComplex=not (args.disabledust), outputName=args.outname, annot=args.annotation, suppressOutputAndFigures=suppressOutputAndFigures)
            # EmAlgo(readsTable, allReadsNum, thresholdTpm=args.tpm, outputName=args.outname, printResult=args.printem, suppressOutputAndFigures=suppressOutputAndFigures)
            print(f"Completing iterations {numTiers - i - 1} tier(s) early because reference tree has {referenceTree.get_num_genes()} genes and was passed {referenceTree.get_num_filters()} filters.",flush=True)

            print("FINAL RESULTS", flush=True)
            print(output_results(emResults, i + 1, capture_gene_group=geneGroup, id_to_old=id_to_old, include_header=False))
            break


if __name__ == '__main__':
    main()
