#!/usr/bin/env python
import re

from src.ManipulateFiles import clean_output, filter_fasta, score_output
from src.CreateMappedReadTable import mapReads
from src.EMstep import EmAlgo
from whichcraft import which
import subprocess as subp
import argparse as argp
import shutil
import sys
import os

__version__ = "1.0.2"

# TODO remove
os.environ['PATH'] = f"/Users/zeliason/Desktop/homebrew/bin:{os.environ.get('PATH')}"

def prereqs():
    programs = ["python3", "samtools", "STAR"]
    ready = True

    for i in range(0, len(programs)):
        if which(programs[i]) is None:
            print(programs[i] + " not installed. Please install " + programs[i])
            ready = False
    return ready



def cmd(args, write=False, filepath=None, verbose=True):
    if(write==True):
        temp = sys.stdout
        sys.stdout = open(filepath, 'w')

        try:
            subp.check_call(args, stdout=subp.DEVNULL if verbose else sys.stdout)
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
            subp.check_call(args, stdout=subp.DEVNULL if verbose else sys.stdout)
        except subp.CalledProcessError as e:
            print("Subprocess error with code: " + str(e.returncode))
            sys.exit(e.returncode)
        except:
            print("An unknown error occurred")
            exit(1)

    return


def filterReferenceFasta(genomeFastaFiles):
    directory, filename = os.path.split(genomeFastaFiles)
    base_name, extension = os.path.splitext(filename)
    if extension.lower() in ('.fa', '.fasta'):
        filteredGenomeFastaFiles = os.path.join(directory, f"{base_name}_ABC{extension}")
    else:
        print(f"The file '{filename}' does not have a .fa or .fasta extension.")

    if not os.path.exists(filteredGenomeFastaFiles):
        filter_fasta(input_file=genomeFastaFiles, output_file=filteredGenomeFastaFiles)

    genomeDir = filteredGenomeFastaFiles + "_STAR"

    return genomeDir, filteredGenomeFastaFiles


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

    # TODO remove shortcut
    myparse.add_argument('--shortcut', action='store_true', default=False)

    # other required arguments
    requiredNamed = myparse.add_argument_group('required arguments')
    requiredNamed.add_argument('-s','--stargenome', help="path to a directory containing STAR-generated\nhuman genome indexes", required=True)

    args = myparse.parse_args()

    if len(args.outname.split('/')) == 1:
        args.outname = os.path.join(os.getcwd(), args.outname)

    reads_dir, reads_basename = os.path.split(args.reads1)
    reads_dir += "/out"

    if not os.path.isdir(args.outname):
        os.makedirs(args.outname)
    base_outname = os.path.basename(args.outname)
    outname = os.path.join(args.outname, base_outname)

    args.starHLA, args.reference = filterReferenceFasta(genomeFastaFiles=args.reference)

    if not os.path.isdir(args.starHLA):
        indexReferenceGenes(genomeDir=args.starHLA, genomeFastaFiles=args.reference, genomeSAindexNbases=args.genomeSAindexNbases, outname=args.outname)

    if args.threads < 1:
        args.threads = 1

    if not prereqs():
        sys.exit(1)


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

    print(args)

    if args.reads2 == "not supplied":
        # TODO remove shortcut
        if not args.shortcut or not os.path.exists(reads_dir + '.Unmapped.out.mate1.fq'):
            print("Aligning reads to human genome")
            cmd(argsHumanAlignSTAR + ["--readFilesIn {}".format(args.reads1)])

            os.rename(reads_dir + '.Unmapped.out.mate1', reads_dir + '.Unmapped.out.mate1.fq')
            os.rename(reads_dir + '.Log.final.out', outname + ".Log.final.out")
        else:
            print('taking shortcut\n')


        with open('{}.Log.final.out'.format(outname),'r') as logFile:
            total_pattern = r"Number of input reads \|	(\d+)"
            unmapped_pattern = r"Number of reads unmapped: ((too short)|(other)) \|\s+(\d+)"

            logFileContents = logFile.read()

            try:
                total_reads = int(re.findall(total_pattern, logFileContents)[0])
                unmapped_reads = sum([int(x[3]) for x in re.findall(unmapped_pattern, logFileContents)])
            except:
                print("ERROR!!!!")
                raise RuntimeError('allReadsNum NOT FOUND!!!!!!!!!!')

            print(f"{unmapped_reads} reads unmapped out of {total_reads} reads total.")

            allReadsNum = total_reads

        if not args.shortcut or not os.path.exists('{}.1.Aligned.out.bam'.format(outname)):
            print("Aligning reads to HLA genomes")
            cmd(["STAR",
                 "--genomeDir {path}".format(path=args.starHLA),
                 "--readFilesIn {sampleName}.Unmapped.out.mate1.fq".format(sampleName=reads_dir),
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
        else:
            print("took second shortcut")

        hlaBams.append('{}.1.Aligned.out.bam'.format(outname))

    else:
        if not args.shortcut or not os.path.exists(reads_dir + '.Unmapped.out.mate1.fq'):
            print("Aligning reads to human genome")
            cmd(argsHumanAlignSTAR + ["--readFilesIn {} {}".format(args.reads1, args.reads2)])

            os.rename(reads_dir + '.Unmapped.out.mate1', reads_dir + '.Unmapped.out.mate1.fq')
            os.rename(reads_dir + '.Unmapped.out.mate2', reads_dir + '.Unmapped.out.mate2.fq')

            os.rename(reads_dir + '.Log.final.out', outname + ".Log.final.out")
        else:
            print('taking shortcut')

        with open('{}.Log.final.out'.format(outname),'r') as logFile:
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
             "--readFilesIn {sampleName}.Unmapped.out.mate1.fq".format(sampleName=reads_dir),
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
        cmd(["STAR", 
             "--genomeDir {path}".format(path=args.starHLA),
             "--readFilesIn {sampleName}.Unmapped.out.mate2.fq".format(sampleName=reads_dir),
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

             "--outFileNamePrefix {}.2.".format(outname)])
        hlaBams.append('{}.2.Aligned.out.bam'.format(outname))

    # Clean unneeded intermediary output and files
    if not args.keepint:
        allowed_extensions = {'.bam', '.tsv', '.pdf', '.csv', '.fq'}

        # Clean both samples and output folders
        directories_to_clean = [os.path.split(reads_dir)[0], args.outname]
        for dir_to_clean in directories_to_clean:
            for filename in os.listdir(dir_to_clean):
                filepath = os.path.join(dir_to_clean, filename)
                if os.path.isdir(filepath):
                    shutil.rmtree(filepath)
                if os.path.isfile(filepath):
                    _, extension = os.path.splitext(filename)
                    if extension.lower() not in allowed_extensions and "Log.final.out" not in filename:
                        os.remove(os.path.join(dir_to_clean, filename))

    print("Creating read table", flush=True)
    readsTable = mapReads(hlaBams, hlaRefPath=args.reference, filterLowComplex=not(args.disabledust), outputName=outname, annot=args.annotation)

    # save readsTable to file
    with open(outname + ".readsTable.tsv", "w") as f:
        for line in readsTable:
            f.write(line + "\n")

    print("Running EM algorithm", flush=True)
    EmAlgo(readsTable, allReadsNum, thresholdTpm=args.tpm, outputName=outname, printResult=args.printem)
    clean_output(outname + ".results.tsv")

    sys.exit(0)


if __name__ == '__main__':
    main()