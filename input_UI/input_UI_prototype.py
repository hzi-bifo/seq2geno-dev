'''
Tasks:
- Parse the arguments
- Determine the working folder
- Check the input files
- Check the dependencies
'''

import snakemake
import os
import argparse

parser = argparse.ArgumentParser(
        description='Seq2Geno: the pipline tool '
            'for genomic features computation')

## functions
parser.add_argument('-v', action= 'version', 
    version='v.Beta')
parser.add_argument('-w', dest='wd', type= str, default= '.',
    help='set the working directory, where the relative paths are counted')
parser.add_argument('-c', dest='comp', action= 'store_true',
    help='compress binary features by pattern')
parser.add_argument('-d', action= 'store_true',
    help='check software dependencies')

## reference
parser.add_argument('--ref_fa', dest='ref_fa', type= str,
    help='reference genome sequences (fasta)', default= '', required=True)
parser.add_argument('--ref_gbk', dest='ref_gbk', type= str,
    help='reference genome annotation (genbank)', default= '', required=True)

## samples
parser.add_argument('--s', dest='samples', type= str,
    help='list of samples and the corresponding reads', default= '', required=True)

## outputs
parser.add_argument('--tree', nargs= '?', dest='tr',
        default= '',  type= str, help='the output tree file')
parser.add_argument('--gpa', nargs= '?', dest='gpa',
        default= '',  type= str, help='the output gene pres/abs table')
parser.add_argument('--s-snp', nargs= '?', dest='s-snp',
        default= '',  type= str, help='the output syn SNPs table')
parser.add_argument('--ns-snp', nargs= '?', dest='ns-snp',
        default= '',  type= str, help='the output non-syn SNPs table')
parser.add_argument('--expr', nargs= '?', dest='xpr',
        default= '',  type= str, help='the output expression table')
parser.add_argument('--ind', nargs= '?', dest='ind',
        default= '',  type= str, help='the output tree table')

args = parser.parse_args()

print(**args)

'''
## create config file
config_txt=''

samples: {}
reference_sequence: {}
reference_annotation: {} 
tree: {}
gpa_table: {}
syn_snps_table: {}
nonsyn_snps_table: {} 
expr_table: {}
indel_table: {}'.format()
'''
'''
w_dir=os.path.abspath('.')
print('working path= {}'.format(w_dir))
seq2geno_dir=os.path.dirname(os.path.realpath(__file__))
print('bin path= {}'.format(seq2geno_dir))
print('config file= {}'.format(args.conf_f))
snakemake.snakemake(snakefile=os.path.join(seq2geno_dir, 'lib/main.smk'),
    configfile= (None if args.conf_f != '' else args.conf_f),
    dryrun=True, workdir= w_dir)


species: 'Paeruginosa'
working_dir: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6'

###### external executives
stampy_exe: '/home/thkuo/bin/stampy-1.0.23/stampy.py'
#raxml_exe: '/home/thkuo/miniconda3/envs/phypal/bin/raxmlHPC-PTHREADS'
#raxml_exe: '/home/thkuo/miniconda3/envs/phypal/bin/raxmlHPC-PTHREADS-AVX2'
raxml_exe: 'raxmlHPC-PTHREADS-AVX2'
#art2genecount_exe: '/net/metagenomics/data/from_moni/old.tzuhao/Paeru/bin/from_SusanneLab/RNA-seq/art2genecount.pl'
#sam2art_exe: '/net/metagenomics/data/from_moni/old.tzuhao/Paeru/bin/from_SusanneLab/RNA-seq/sam2art.pl'

software:
  mapper: 'bwa'
  assembler: 'spades'

cores: 16

##### files
#### sample wise
samples: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/data/samples.tsv'
#### only one
reference_sequence: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/data/reference/RefCln_UCBPP-PA14.fa'
reference_annotation: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/data/reference/RefCln_UCBPP-PA14.gbk'
##roary_gpa:
tree: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/results/Paeru.nwk'
gpa_table: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/results/gpa.tab'
syn_snps_table: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/results/syn_snps.tab'
nonsyn_snps_table: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/results/nonsyn_snps.tab'
expr_table: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/results/expr.tab'
indel_table: '/net/metagenomics/data/from_moni/old.tzuhao/seq2geno/dev_versions/v6/results/annot.tab'
'''
