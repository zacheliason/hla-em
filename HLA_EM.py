from src.ManipulateFiles import predict_genotype_from_MLE, filter_fasta, plot_coverage_maps, plot_pie_charts
from src.CreateMappedReadTable import mapReads
from src.EMstep import EmAlgo
from whichcraft import which
import subprocess as subp
import argparse as argp
import traceback
import shutil
import time
import sys
import os
import re


__version__ = "1.0.2"

# TODO remove
os.environ['PATH'] = f"/Users/zeliason/Desktop/homebrew/bin:{os.environ.get('PATH')}"

def prereqs():
    programs = ["python3", "samtools"]#, "STAR"]
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

def get_read_counts_from_log(filepath):
    with open(filepath, 'r') as logFile:
        total_pattern = r"Number of input reads \|	(\d+)"
        unmapped_pattern = r"Number of reads unmapped: ((too short)|(other)) \|\s+(\d+)"

        log_contents = logFile.read()

        try:
            total_reads = int(re.findall(total_pattern, log_contents)[0])
            unmapped_reads = sum([int(x[3]) for x in re.findall(unmapped_pattern, log_contents)])
        except:
            raise RuntimeError('allReadsNum not found')

        return total_reads, unmapped_reads


def main(args=None):
    if args is None:
        installDir = os.path.dirname(os.path.abspath(__file__))

        myparse = argp.ArgumentParser(prog='HLA-EM', description='HLA-EM is an HLA genotyping tool that utilizes an expectation maximization algorithm to identify the presence of different HLA genotypes in a sample from RNA-seq data.', formatter_class=lambda prog: argp.RawTextHelpFormatter(prog, width=99999))

        # positional arguments
        myparse.add_argument("reads1", help="single-end FASTQ file or first paired-end FASTQ file")
        myparse.add_argument("reads2", nargs='?', help="(optional) second paired-end FASTQ file", default="not supplied")

        # options
        myparse.add_argument('-t','--threads', type=int,  help="number of threads to use [1]", default=1)
        myparse.add_argument('-g','--genomeSAindexNbases', type=int,  help="number of bases to use [6]", default=6)
        myparse.add_argument('-r','--reference', help="HLA reference genome in FASTA format,\nto be used in place of default HLA reference", default='hla_gen.fasta')
        myparse.add_argument('-a','--annotation', help="HLA gene annotations in TSV format,\nto be used in place of default HLA annotations\n[{}]".format('$HLA-EMPath/reference/hla_gene_annot.tsv'), default=installDir+'/reference/hla_gene_annot.tsv')
        myparse.add_argument('--starHLA', help="path to a directory containing STAR-generated\nHLA genome indexes based on the above FASTA", default=0)
        myparse.add_argument('-o', '--outname', type=str, help="output file name prefix [./hlaEM]", default='./hlaEM')
        myparse.add_argument('-d', '--disabledust', action='store_true', help="disable filtering of low-complexity reads")
        myparse.add_argument('--tpm', type=float, help="TPM threshold for identifying a true positive [1.48]", default=1.48)
        myparse.add_argument('-p', '--printem', action='store_true', help="print EM results to STDOUT")
        myparse.add_argument('-k', '--keepint', action='store_true', help="keep intermediate files", default=False)
        myparse.add_argument('--suppress_figs', action='store_true', help="skip plots for faster performance")
        myparse.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))
        myparse.add_argument('--training', type=str, default="")


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

    # TODO remove shortcut
    if args.shortcut and os.path.exists(os.path.join(args.outname, 'final_predictions.csv')):
        pass
        # print(f"final predictions already exists! Taking shortcut :)")
        # exit(0)

    if not os.path.isdir(args.outname):
        os.makedirs(args.outname)
    base_outname = os.path.basename(args.outname)
    outname = os.path.join(args.outname, base_outname)

    args.starHLA, args.reference = filterReferenceFasta(genomeFastaFiles=args.reference)
    if not os.path.exists(args.starHLA):
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

    if args.reads2 == "not supplied":
        # TODO remove shortcut
        if not args.shortcut or not (os.path.exists(reads_dir + '.Unmapped.out.mate1.fastq') and os.path.exists('{}.Log.final.out'.format(outname))):
            print("Aligning reads to human genome")
            cmd(argsHumanAlignSTAR + ["--readFilesIn {}".format(args.reads1)])

            os.rename(reads_dir + '.Unmapped.out.mate1', reads_dir + '.Unmapped.out.mate1.fastq')
            os.rename(reads_dir + '.Log.final.out', outname + ".Log.final.out")

        allReadsNum, unmappedReadsNum = get_read_counts_from_log(outname + ".Log.final.out")
        print(f"{unmappedReadsNum} reads unmapped out of {allReadsNum} reads total")

        # TODO remove shortcut
        if not args.shortcut or not os.path.exists('{}.1.Aligned.out.bam'.format(outname)):
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
        # TODO remove shortcut
        if not args.shortcut or not (os.path.exists(reads_dir + '.Unmapped.out.mate1.fastq') and os.path.exists('{}.Log.final.out'.format(outname)) and os.path.exists(os.path.exists(reads_dir + '.Unmapped.out.mate2.fastq'))):
            print("Aligning paired reads to human genome")
            cmd(argsHumanAlignSTAR + ["--readFilesIn {} {}".format(args.reads1, args.reads2)])

            os.rename(reads_dir + '.Unmapped.out.mate1', reads_dir + '.Unmapped.out.mate1.fastq')
            os.rename(reads_dir + '.Unmapped.out.mate2', reads_dir + '.Unmapped.out.mate2.fastq')

            os.rename(reads_dir + '.Log.final.out', outname + ".Log.final.out")

        allReadsNum, unmappedReadsNum = get_read_counts_from_log(outname + ".Log.final.out")
        # print(f"{unmappedReadsNum} reads unmapped out of {allReadsNum} reads total")

        # TODO remove shortcut
        if not args.shortcut or not (os.path.exists('{}.1.Aligned.out.bam'.format(outname))):
            print("Aligning reads to HLA genomes")

            cmd(["STAR",
                 "--genomeDir {path}".format(path=args.starHLA),
                 "--readFilesIn {sampleName}.Unmapped.out.mate1.fastq {sampleName}.Unmapped.out.mate2.fastq".format(
                     sampleName=reads_dir),
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

    # Clean unneeded intermediary output and files
    if not args.keepint:
        allowed_extensions = {'.bam', '.tsv', '.pdf', '.csv', '.fastq', '.fq', '.json'}

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

    if not args.shortcut or not (os.path.exists('{}.mappedReads.tsv'.format(outname))):
        print("Creating read table", flush=True)
        time_start = time.time()
        readsTable = mapReads(hlaBams, hlaRefPath=args.reference, filterLowComplex=not(args.disabledust), outputName=outname, annot=args.annotation, suppressOutputAndFigures=args.suppress_figs)
        time_end = time.time()

        print(f"CreateMappedReads took {time_end - time_start} seconds")

    if not args.shortcut or not (os.path.exists('{}.results.tsv'.format(outname))):
        print("Running EM algorithm", flush=True)
        if args.shortcut:
            with open('{}.mappedReads.tsv'.format(outname)) as f:
                readsTable = f.read().split('\n')

        time_start = time.time()
        EmAlgo(readsTable, outname=outname, thresholdTpm=args.tpm)
        time_end = time.time()
        print(f"EM algorithm took {time_end - time_start} seconds")

    predictions = predict_genotype_from_MLE(outname + ".results.tsv", outname, base_outname, training_spreadsheet=args.training)
    predicted_types = predictions.index.values.tolist()

    print("PREDICTED HLA TYPES:")
    for predicted_type in predicted_types:
        print(" " + predicted_type)
    print()

    if not args.suppress_figs:
        # TODO remove try/except
        try:
            plot_pie_charts(outname + ".final_predictions.csv", outname + ".results.tsv", outname)
            # plot_coverage_maps(outname + ".cov_plot_args.json", predicted_types)
        except:
            print(traceback.format_exc())
            pass

    # TODO restore
    # sys.exit(0)


if __name__ == '__main__':
    main()
