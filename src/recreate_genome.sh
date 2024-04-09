#!/bin/bash

# Based on instructions from @alexdobin
# https://github.com/alexdobin/STAR/issues/317

# Transform characters and output to G1 file
cat Genome | tr '\0\1\2\3\4' 'ACGTN' | tr -s '\5' '\n' > G1

# Count lines in G1 and chrName.txt (G1 may contains many more lines if you included GTF or splice junctions in the genome generation.)
wc -l G1 chrName.txt

# Extract the number of lines in chrName.txt
n=$(wc -l chrName.txt | awk '{print $1}')

# Output the first n lines of G1 into G2
head -n $n G1 > G2

# Prepend '>' to each line in chrName.txt and output to chN
sed 's/^/>/' chrName.txt > chN

# Combine chN and G2 line by line and output to Genome.fa
paste -d'\n' chN G2 > Genome.fa

# Clean up
rm G1 G2 chN