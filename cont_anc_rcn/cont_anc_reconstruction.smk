#TREE_OUT=config['TREE_OUT']
#EXPR_OUT=config['EXPR_OUT']
import os
import yaml
expr_config_f= config['expr_config']
phylo_config_f= config['phylo_config']
expr_config= yaml.load(open(expr_config_f, 'r'))
phylo_config= yaml.load(open(phylo_config_f, 'r'))
EXPR_OUT= os.path.join(os.path.dirname(expr_config_f), expr_config['out_table'])
TREE_OUT= os.path.join(os.path.dirname(phylo_config_f), phylo_config['tree_f'])
C_ANCREC_OUT= config['C_ANCREC_OUT']
rule all:
    input:
        output_dir=C_ANCREC_OUT
        
rule cont_anc_reconstruction:
    input:
        tree_f= TREE_OUT,
        data_f= EXPR_OUT
    output: 
        output_dir=directory(C_ANCREC_OUT)
    conda: 'ar_env.yaml'
    script: 'contAncRec.R'
