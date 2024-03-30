import HLA_EM
import argparse
import os

# Define the volume path
VOLUME_PATH = "/Users/zacheliason/Downloads/hla-em"  # Adjust this if needed

# Define the range of trials (0 to 31)
trials = range(32)

# Function to simulate argparse and create namespace
def create_namespace():
    parser = argparse.ArgumentParser()

    parser.add_argument("reads1", help="single-end FASTQ file or first paired-end FASTQ file")
    parser.add_argument("reads2", nargs='?', help="(optional) second paired-end FASTQ file", default="not supplied")

    parser.add_argument('-t', '--threads', type=int, help="number of threads to use [1]", default=1)
    parser.add_argument('-g', '--genomeSAindexNbases', type=int, help="number of bases to use [6]", default=6)
    parser.add_argument('-r', '--reference',
                         help="HLA reference genome in FASTA format,\nto be used in place of default HLA reference",
                         default='hla_gen.fasta')
    parser.add_argument('--starHLA',
                         help="path to a directory containing STAR-generated\nHLA genome indexes based on the above FASTA",
                         default=0)
    parser.add_argument('-o', '--outname', type=str, help="output file name prefix [./hlaEM]", default='./hlaEM')
    parser.add_argument('--annotation', type=str, help="output file name prefix [./hlaEM]", default='./hlaEM')
    parser.add_argument('-d', '--disabledust', action='store_true', help="disable filtering of low-complexity reads")
    parser.add_argument('--tpm', type=float, help="TPM threshold for identifying a true positive [1.48]", default=1.48)
    parser.add_argument('-p', '--printem', action='store_true', help="print EM results to STDOUT")
    parser.add_argument('-k', '--keepint', action='store_true', help="keep intermediate files", default=False)
    parser.add_argument('--suppress_figs', action='store_true', help="skip coverage plots for faster performance")
    parser.add_argument('--training', type=str, default="")
    parser.add_argument('--shortcut', action='store_true', default=False)

    # other required arguments
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-s', '--stargenome',
                               help="path to a directory containing STAR-generated\nhuman genome indexes",
                               required=True)

    return parser.parse_args([
        "-t", "4",
        "--training", f"{OUT_DIR}/training.tsv",
        "--shortcut",
        "-o", OUTPUT_DIR,
        "-s", os.path.join(VOLUME_PATH, "EnsembleGenome_STAR_without_scaffolds"),
        "-r", os.path.join(VOLUME_PATH, "hla_gen.fasta"),
        os.path.join(TRIAL_DIR, "sim.HLA.reads_01.fq"),
        os.path.join(TRIAL_DIR, "sim.HLA.reads_02.fq")
    ])

# Iterate over each trial
for trial in trials:
    TRIAL_DIR = os.path.join(VOLUME_PATH, f"reference/samples/trial_{trial}")
    OUT_DIR = os.path.join(VOLUME_PATH, f"output_paired")
    OUTPUT_DIR = os.path.join(OUT_DIR, f"trial_{trial}")

    # Create namespace with arguments
    args = create_namespace()

    # Call HLA_EM.main with arguments
    HLA_EM.main(args)