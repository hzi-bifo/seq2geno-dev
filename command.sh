#!/bin/bash
#$ -V
#$ -l h_core=20
#$ -l h_vmem=100G

cd /net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6
source activate seq2geno

#snakemake --unlock --snakefile=MAIN.smk
#snakemake --notemp --snakefile=MAIN.smk
#snakemake --notemp -R load_reference --snakefile=test.smk
#snakemake --notemp --snakefile=test.smk
#snakemake --notemp -R rna_mapping --snakefile=MAIN.smk seq2geno_temp/CH2500/rna.sin
snakemake --notemp -R load_reference --snakefile=MAIN.smk
