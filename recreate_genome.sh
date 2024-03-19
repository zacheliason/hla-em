cat Genome | tr '\0\1\2\3\4' 'ACGTN' | tr -s '\5' '\n' > G1
wc -l G1 chrName.txt
n=wc -l chrName.txt | cut -d' ' -f1
head -n $n G1 > G2
sed 's/^/>/' chrName.txt > chN
paste -d'\n' chN G2 > Genome.fa
