from CreateMappedReadTable import mapReads

import subprocess
import os

# Define the SAMtools command to execute
os.environ['PATH'] = '/Users/zeliason/Desktop/homebrew/bin:' + os.environ['PATH']
os.system('echo $PATH')



# # Run the command and capture the output
# try:
#     output = subprocess.check_output(command, stderr=subprocess.STDOUT, text=True)
#     print("SAMtools is installed and working. Here's the version information:")
#     print(output)
# except subprocess.CalledProcessError as e:
#     print("An error occurred while running SAMtools command:")
#     print(e.output)
#
command = ['samtools', '--version']

# Run the command and capture the output
try:
    output = subprocess.check_output(command, stderr=subprocess.STDOUT, text=True)
    print("SAMtools is installed and working. Here's the version information:")
    print(output)
except subprocess.CalledProcessError as e:
    print("An error occurred while running SAMtools command:")
    print(e.output)

# def wmapReads(hlaBams, defaultHlaRef=True, hlaRefPath='', annot='', filterLowComplex=True, outputName='hlaType', covMapYmax=0, suppressOutputAndFigures=True):
look = mapReads(hlaBams=['/Users/zeliason/Desktop/hla-em/src/output/trial_0/trial_0.1.Aligned.out.bam'], hlaRefPath='/Users/zeliason/Desktop/hla-em/src/hla_gen.fasta', filterLowComplex=True)