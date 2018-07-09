'''
All the files listed in this interface should be specified by the user (in config.yaml or samples.tsv),
no matter a file exists or not.
If a file is already there, let the following workflow use it. 
Otherwise, compute that file with the tools chosen by the user.

- List all the files required, allow the users to specify them, and this tool should be able to parse check them => follow the three execution phases of snakemake
- This tools cares the required files and should be influenced care as few as possible about the software otpions and usages => a genome written by spades.py or by hand should both be acceptable by the following workflow
- Should the detailed workflows (eg. the exact running parts) by included (stronger dependency) or sub-workflow (another working environment and config file required)
'''

import pandas as pd
import os

configfile: "config.yaml"
samples_df=pd.read_table(config["samples"], sep= '\t', header= 0).set_index("strain", drop=False)
strains=samples_df['strain'].tolist()

rule all:
    input:
        ml_markdown= config['ml_markdown']

rule import_dna_reads:
    input:  
    output:
        dna_reads

rule import_dna_params:
    input:

    output:
        dna_params

rule import_reference:
    input:

    output:
        ref

rule import_rna_reads:
    input:
    
    output:
        rna_reads

rule import_rna_params:
    input:

    output:
        rna_params

rule count_gpa:
    input:
        dna_reads
        dna_params
    output:
        roary_gpa

rule detectSNPs:
    input:
        dna_reads
        dna_params
        ref
    output:
        vcf

rule compute_expressions:
    input:
        rna_reads
        rna_params
        ref
    output:
        

rule create_gpa_table:
    input:
        roary_gpa
    output:
        gpa_table

rule create_snps_table:
    input:
        vcfs
    output:
        snps_table= config['snps_table']

rule infer_tree:
    input:
        vcfs
        ref
    output:
        tree= config['tree']

rule create_expr_table:
    input:
        
    output:
        expr_table= config['expr_table']
             
rule make_ml_input:
    input:
        gpa_table= config['gpa_table']
        snps_table= config['snps_table']
        tree= config['tree']
        expr_table= config['expr_table']
    output:
        ml_markdown= config['ml_markdown']
