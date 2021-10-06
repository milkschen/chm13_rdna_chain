set -x

# Required files:
#  - chm13.draft_v1.0.fasta
#      CHM13 v1.0 reference
#      from https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.v1.0.fasta
#      and  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.v1.0.fasta.fai
#  - chm13.v1.1.fasta
#      CHM13 v1.1 reference
#      from https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.v1.1.fasta
#      and  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.v1.1.fasta.fai
#  - v1.0.list
#      rDNA flanking sequence regions in CHM13 v1.0
#      from flanking/v1.0.list
#  - v1.1.list
#      rDNA flanking sequence regions in CHM13 v1.1
#      from flanking/v1.1.list

wget -nc https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.v1.0.fasta
wget -nc https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.v1.0.fasta.fai
wget -nc https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.v1.1.fasta
wget -nc https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.v1.1.fasta.fai

# Extract flanking sequences
samtools faidx -r flanking/v1.0.list chm13.v1.0.fasta > chm13.v1.0.rDNA_flanking.fa
samtools faidx -r flanking/v1.1.list chm13.v1.1.fasta > chm13.v1.1.rDNA_flanking.fa

# Align flanking sequences
meryl count k=19 output merylDB_v1.0 chm13.v1.0.fasta
meryl print greater-than distinct=0.9998 merylDB_v1.0 > repetitive_k19_v1.0.txt
winnowmap -k 19 -c -W repetitive_k19_v1.0.txt chm13.v1.0.fasta chm13.v1.1.rDNA_flanking.fa > chm13.v1.1.rDNA_flanking.to_v1.0.paf
winnowmap -k 19 -a -W repetitive_k19_v1.0.txt chm13.v1.0.fasta chm13.v1.0.rDNA_flanking.fa > chm13.v1.0.rDNA_flanking.to_v1.0.sam

meryl count k=19 output merylDB_v1.1 chm13.v1.1.fasta
meryl print greater-than distinct=0.9998 merylDB_v1.1 > repetitive_k19_v1.1.txt
winnowmap -k 19 -a -W repetitive_k19_v1.1.txt chm13.v1.1.fasta chm13.v1.1.rDNA_flanking.fa > chm13.v1.1.rDNA_flanking.to_v1.1.sam
winnowmap -k 19 -c -W repetitive_k19_v1.1.txt chm13.v1.1.fasta chm13.v1.0.rDNA_flanking.fa > chm13.v1.0.rDNA_flanking.to_v1.1.paf

# Build chain from 1.1 to 1.0
wget -nc https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/changes/v1.0_to_v1.1/v1.1_to_v1.0.chain
python rdna_paf_to_chain.py -p chm13.v1.1.rDNA_flanking.to_v1.0.paf  -a chm13.v1.1.rDNA_flanking.to_v1.1.sam  -o chm13.v1.1.rDNA_flanking.to_v1.0.chain -ic 25
cp v1.1_to_v1.0.chain v1.1_to_v1.0_rdna.chain
cat chm13.v1.1.rDNA_flanking.to_v1.0.chain >> v1.1_to_v1.0_rdna.chain
python merge_colinear_chain.py -c v1.1_to_v1.0_rdna.chain -o v1.1_to_v1.0_rdna.merged.chain

# Build chain from 1.0 to 1.1
wget -nc https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/changes/v1.0_to_v1.1/v1.0_to_v1.1.chain
python rdna_paf_to_chain.py -p chm13.v1.0.rDNA_flanking.to_v1.1.paf -a chm13.v1.0.rDNA_flanking.to_v1.0.sam -o chm13.v1.0.rDNA_flanking.to_v1.1.chain -ic 25
cp v1.0_to_v1.1.chain v1.0_to_v1.1_rdna.chain
cat chm13.v1.0.rDNA_flanking.to_v1.1.chain >> v1.0_to_v1.1_rdna.chain
python merge_colinear_chain.py -c v1.0_to_v1.1_rdna.chain -o v1.0_to_v1.1_rdna.merged.chain
