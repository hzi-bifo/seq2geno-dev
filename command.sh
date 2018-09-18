#!/bin/bash
#$ -V
#$ -l h_core=30
#$ -l h_vmem=100G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6
source activate seq2geno

#snakemake --unlock --snakefile=MAIN.smk
#snakemake --notemp --snakefile=MAIN.smk
snakemake --notemp -R rna_mapping --snakefile=test.smk
