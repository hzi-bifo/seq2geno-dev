#!/bin/bash
#$ -V
#$ -l h_core=30
#$ -l h_vmem=120G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6
source activate seq2geno
source activate Ariane_dna

#/home/thkuo/bin/stampy-1.0.23/stampy.py --bwaoptions="-q10 data/reference/RefCln_UCBPP-PA14.fa" -g data/reference/RefCln_UCBPP-PA14.fa -h data/reference/RefCln_UCBPP-PA14.fa -t30 -M data/strains/rna/CH2500.fastq.gz > ./seq2geno_temp/CH2500/rna.sam
/home/thkuo/bin/stampy-1.0.23/stampy.py --bwaoptions="-q10 /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/data/reference/RefCln_UCBPP-PA14.fa"         -g /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/data/reference/RefCln_UCBPP-PA14.fa -h /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/data/reference/RefCln_UCBPP-PA14.fa -t30         -M data/strains/rna/CH2500.fastq.gz > seq2geno_temp/CH2500/stampy/rna.sam
