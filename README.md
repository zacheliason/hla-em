# HLA-EM

HLA-EM is an HLA genotyping tool that utilizes an expectation maximization algorithm to identify a patient's HLA genotype in a sample from RNA-seq data. 

## Prerequisites
  - [python](https://www.python.org/) (v 2.7 or 3+)
    - [NumPy](http://http://www.numpy.org/)
    - [Matplotlib](https://matplotlib.org/)
    - [whichcraft](https://pypi.org/project/whichcraft/)
  - [STAR](https://github.com/alexdobin/STAR)
    - to compile STAR, you may need to install [CMAKE](https://cmake.org/) or [gcc](https://gcc.gnu.org/)
  - [SAMtools](http://samtools.sourceforge.net/)

Alternatively, you can use a Docker image that includes all the necessary dependencies. The Docker image eliminates the need for manual installations and ensures a consistent environment. To use the Docker image, follow the steps below:

### Using Docker Image
  1. Install [Docker](https://www.docker.com/) on your machine.
  2. Pull the Docker image:

     ```bash
     docker pull zeliason/hla_em:latest
     ```

  3. Run your application within the Docker container:

     ```bash
     docker run -it --rm -v /path/to/hla-em/src:/src zeliason/hla_em:latest /bin/bash
     ```

  4. You may now run the HLA_EM.py as normal within the container from the mounted `/src` volume.
  
## Installation
First ensure the prerequisites have been installed on your system and that STAR and SAMtools appear in your PATH.  In order for STAR to align reads against the human genome, you will first need to obtain the human reference genome in FASTA format and a corresponding annotation file in GTF format (these can be downloaded from [Ensembl](https://ensembl.org/Homo_sapiens/Info/Index)).  You must then generate genome indexes for STAR aligner using its genomeGenerate command (see the STAR documentation for details).

You can then download HLA-EM from https://github.com/zacheliason/hla-em.  Click on "Clone or download" and then click on "Download ZIP".  Unzip the contents in the desired location.  
If you desire, you can include HLA-EM in your PATH by running:
```
$ export PATH=$INSTALLDIR:$PATH
```
where `$INSTALLDIR` is the folder into which you extracted the HLA-EM files.
  
  
## Input
Running HLA-EM with the -h option (or --help) will print a desciption of its optional and required input arguments.  A description of each follows.
```
HLA-EM.py [-h] [-t THREADS] [-r REFERENCE] [--starHLA STARHLA] [-o OUTNAME]
          [-d] [--tpm TPM] [-p] [-k] [-v] -s STARGENOME reads1 [reads2]
```
### Optional arguments
- -h, --help  
      Prints the help menu, which contains an example of the function usage and abbreviated explanations of each of the options.
- -t THREADS, --threads THREADS  
     Some portions of the HLA-EM pipeline support multithreading; if you wish to use this feature, set this to the number of cores available (default: 1)
- -r REFERENCE, --reference REFERENCE  
     HLA-EM includes a comprehensive set of HLA reference genes maintained by [IMGT, the international ImMunoGeneTics Database](https://github.com/ANHIG/IMGTHLA/blob/Latest/hla_gen.fasta), the use of which is recommended.  However, if you wish to supply your own HLA gene reference file in FASTA format, use this option with the path to the file
- -a ANNOTATION, --annotation ANNOTATION  
     A file of HLA gene annotations in TSV format to be used in place of the default HLA gene annotations [$INSTALLDIR/reference/hla_gene_annot.tsv]
- --starHLA STARHLA
     path to a directory containing STAR-generated HLA genome indexes based on the reference FASTA
- -o OUTNAME, --outname OUTNAME  
     Prefix attached to output file names (default: ./hlaEM)
- -d, --disabledust  
     By default, HLA-EM filters low-complexity reads using the DUST algorithm.  Specify this option to disable filtering
- --tpm TPM  
     TPM threshold for identifying a true positive HLA genotype (default: 1.48)
- -p, --printem  
     Causes HLA-EM to print its results to STDOUT
- -k, --keepint  
     By default, HLA-EM removes files generated by intermediate steps of the pipeline. Specify this option to keep these files
- -v, --version  
     Print program's version number and exit

### Required arguments
- -s STARGENOME, --stargenome STARGENOME  
     Path to the directory containing STAR-generated human genome index files

### Positional arguments
- reads1 : A FASTQ file of single-end RNA-seq reads or the first of two paired-end FASTQs 
- reads2 : (optional) Second paired-end FASTQ file  
  Reads files may be gzipped; if so, the filename must end in ".gz".

## Output
  The output of the EM algorithm is written to OUTNAME.results.tsv (and, optionally, printed to stdout).  The first two lines indicate how may steps the algorithm took to converge and the maximum likelihood estimate (MLE) of the sequencing error rate.  The remaining lines have the format:
  
```
HLAtype   MappedReads   MappedProportion   MLE_Reads   MLE_Probability
```
This is a table including each HLA type with at least one read mapped to it, indicating the number of mapped reads, the proportion of all reads mapped to the HLA reference genomes mapped to this type, the MLE of the expected number of reads for this HLA type, and the MLE probability of this HLA type.

In addition to the results table, HLA-EM also generates read coverage maps for each HLA type with a non-zero MLE probability (OUTNAME.\*.cov.pdf) and a visualization of the difference between the mapped reads proportions and MLE probabilites (OUTNAME.props.pdf).  The reads aligned to the HLA references are also recorded in BAM format (OUTNAME.aligned.\*.bam).
