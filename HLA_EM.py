#!/usr/bin/env python

from CreateMappedReadTable import mapReads
from EMstep import EmAlgo
from whichcraft import which
import subprocess as subp
import argparse as argp
import shutil
import sys
import os

__version__ = "1.0.2"


def prereqs():
    programs = ["python3", "samtools", "STAR"]
    ready = True

    for i in range(0, len(programs)):
        if which(programs[i]) is None:
            print(programs[i] + " not installed. Please install " + programs[i])
            ready = False
    return ready


def cmd(args, write=False, filepath=None):
    if(write==True):
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
    
    myparse = argp.ArgumentParser(prog='HLA-EM', description='HLA-EM is an HLA genotyping tool that utilizes an expectation maximization algorithm to identify the presence of different HLA genotypes in a sample from RNA-seq data.', formatter_class=lambda prog: argp.RawTextHelpFormatter(prog, width=99999))
    
    # positional arguments
    myparse.add_argument("reads1", help="single-end FASTQ file or first paired-end FASTQ file")
    myparse.add_argument("reads2", nargs='?', help="(optional) second paired-end FASTQ file", default="not supplied")

    # options
    myparse.add_argument('-t','--threads', type=int,  help="number of threads to use [1]", default=1)
    myparse.add_argument('-g','--genomeSAindexNbases', type=int,  help="number of bases to use [6]", default=6)
    myparse.add_argument('-r','--reference', help="HLA reference genome in FASTA format,\nto be used in place of default HLA reference", default=0)
    myparse.add_argument('-a','--annotation', help="HLA gene annotations in TSV format,\nto be used in place of default HLA annotations\n[{}]".format('$HLA-EMPath/reference/hla_gene_annot.tsv'), default=installDir+'/reference/hla_gene_annot.tsv')
    myparse.add_argument('--starHLA', help="path to a directory containing STAR-generated\nHLA genome indexes based on the above FASTA", default=0)
    myparse.add_argument('-o', '--outname', type=str, help="output file name prefix [./hlaEM]", default='./hlaEM')
    myparse.add_argument('-d', '--disabledust', action='store_true', help="disable filtering of low-complexity reads")
    myparse.add_argument('--tpm', type=float, help="TPM threshold for identifying a true positive [1.48]", default=1.48)
    myparse.add_argument('-p', '--printem', action='store_true', help="print EM results to STDOUT")
    myparse.add_argument('-k', '--keepint', action='store_true', help="keep intermediate files")
    myparse.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))

    # other required arguments
    requiredNamed = myparse.add_argument_group('required arguments')
    requiredNamed.add_argument('-s','--stargenome', help="path to a directory containing STAR-generated\nhuman genome indexes", required=True)

    args = myparse.parse_args()

    # TODO remove defaultHlaRef
    defaultHlaRef = False
    if args.starHLA == 0:
        args.starHLA = args.reference + '_STAR'
        if not os.path.isdir(args.starHLA):
            indexReferenceGenes(genomeDir=args.starHLA, genomeFastaFiles=args.reference,
                                genomeSAindexNbases=args.genomeSAindexNbases, outname=args.outname)

    if args.starHLA == 0:
        print('Please provide the path to a folder of STAR indices based on your specified HLA genome using the --starHLA argument')
        sys.exit(1)

    if args.threads < 1:
        args.threads = 1

    if not prereqs():
        sys.exit(1)

    if len(args.outname.split('/')) == 1:
        args.outname = os.path.join(os.getcwd(), args.outname)

    reads_basename, reads_dir = os.path.split(args.reads1)
    print(reads_dir)

    if not os.path.isdir(args.outname):
        cmd(["mkdir", args.outname])
    base_outname = os.path.basename(args.outname)
    outname = os.path.join(args.outname, base_outname)

    allReadsNum = -1
    hlaBams = []
    argsHumanAlignSTAR = ["STAR", 
             "--genomeDir {path}".format(path=args.stargenome),
             "--runThreadN {}".format(args.threads),
             "--chimSegmentMin 18",
             "--outSAMtype BAM Unsorted",
             "--outReadsUnmapped Fastx",
             "--outFilterMultimapNmax 100",
             "--outFilterMismatchNmax 4",
             "--outFileNamePrefix {}.".format(reads_dir)]
    if args.reads1.endswith(".gz"):
        argsHumanAlignSTAR.append("--readFilesCommand zcat")
    if args.reads2 == "not supplied":
        print("Aligning reads to human genome")
        cmd(argsHumanAlignSTAR + ["--readFilesIn {}".format(args.reads1)])

        with open('{}.Log.final.out'.format(reads_dir),'r') as logFile:
            for line in logFile:
                line = line.strip()
                if line.startswith('Number of input reads'):
                    allReadsNum = int(line.split()[-1])
                    print("Total reads: {}".format(allReadsNum))
                    break

        print("Aligning reads to HLA genomes")
        cmd(["STAR", 
             "--genomeDir {path}".format(path=args.starHLA),
             "--readFilesIn {sampleName}.Unmapped.out.mate1.fastq".format(sampleName=reads_dir),
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
             "--outFileNamePrefix {}.1.".format(outname)])
        hlaBams.append('{}.1.Aligned.out.bam'.format(outname))

    else:
        print("Aligning reads to human genome")
        cmd(argsHumanAlignSTAR + ["--readFilesIn {} {}".format(args.reads1, args.reads2)])

        with open('{}.Log.final.out'.format(reads_dir),'r') as logFile:
            for line in logFile:
                line = line.strip()
                if line.startswith('Number of input reads'):
                    allReadsNum = int(line.split()[-1])
                    print(f"allReadsNum: {allReadsNum}")
                    print(line)
                    print()
                    break

        print("Aligning reads to HLA genomes")
        cmd(["STAR", 
             "--genomeDir {path}".format(path=args.starHLA),
             "--readFilesIn {sampleName}/unmapped_reads1.fastq".format(sampleName=reads_dir),
             "--runThreadN {}".format(args.threads),
             "--twopassMode Basic",
             "--outSAMtype BAM Unsorted",
             "--outSAMattributes NH HI NM MD AS XS",
             "--outFilterMultimapNmax 999",
             "--outFilterMismatchNmax 999",
             "--outFilterMismatchNoverLmax 0.08",
             "--outFileNamePrefix {}.1.".format(outname)])
        hlaBams.append('{}.1.Aligned.out.bam'.format(outname))
        cmd(["STAR", 
             "--genomeDir {path}".format(path=args.starHLA),
             "--readFilesIn {sampleName}/unmapped_reads2.fastq".format(sampleName=reads_dir),
             "--runThreadN {}".format(args.threads),
             "--twopassMode Basic",
             "--outSAMtype BAM Unsorted",
             "--outSAMattributes NH HI NM MD AS XS",
             "--outFilterMultimapNmax 999",
             "--outFilterMismatchNmax 999",
             "--outFilterMismatchNoverLmax 0.08",
             "--outFileNamePrefix {}.2.".format(outname)])
        hlaBams.append('{}.2.Aligned.out.bam'.format(outname))

    if not args.keepint:
        os.remove('{}.Log.progress.out'.format(outname))
        os.remove('{}.Log.final.out'.format(outname))
        os.remove('{}.Log.out'.format(outname))
        os.remove('{}.SJ.out.tab'.format(outname))
        os.remove('{}.Chimeric.out.junction'.format(outname))
        os.remove('{}.Aligned.out.bam'.format(outname))
        # os.remove('{}.Unmapped.out.mate1'.format(outname))

        os.remove('{}.1.Log.progress.out'.format(outname))
        os.remove('{}.1.Log.final.out'.format(outname))
        os.remove('{}.1.Log.out'.format(outname))
        os.remove('{}.1.SJ.out.tab'.format(outname))
        shutil.rmtree('{}.1._STARgenome'.format(outname))
        shutil.rmtree('{}.1._STARpass1'.format(outname))
        if args.reads2 != "not supplied":
            # os.remove('{}.Unmapped.out.mate2'.format(outname))
            os.remove('{}.2.Log.progress.out'.format(outname))
            os.remove('{}.2.Log.final.out'.format(outname))
            os.remove('{}.2.Log.out'.format(outname))
            os.remove('{}.2.SJ.out.tab'.format(outname))
            shutil.rmtree('{}.2._STARgenome'.format(outname))
            shutil.rmtree('{}.2._STARpass1'.format(outname))

    print("Creating read table", flush=True)
    readsTable = mapReads(hlaBams, defaultHlaRef=defaultHlaRef, hlaRefPath=args.reference, filterLowComplex=not(args.disabledust), outputName=outname, annot=args.annotation)

    # save readsTable to file
    with open(outname + ".readsTable.tsv", "w") as f:
        for line in readsTable:
            f.write(line + "\n")


    print("Running EM algorithm", flush=True)
    EmAlgo(readsTable, allReadsNum, thresholdTpm=args.tpm, outputName=outname, printResult=args.printem)

    sys.exit(0)


if __name__ == '__main__':
    main()