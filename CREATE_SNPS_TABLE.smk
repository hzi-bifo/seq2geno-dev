#####
### Ariane's pipeline
### Remember to include all the requird scripts
### Problems:
### 1. where and how to include the required scripts?
### 2. python2
### 3. requiring flatcounts


rule covert_snps_table:
    input:
        non_syn_output= os.path.join(TMP_D, 'non-syn_SNPs.tab'),
        syn_output= os.path.join(TMP_D, 'syn_SNPs.tab')
    output: 
        config['nonsyn_snps_table'],
        config['syn_snps_table']
    params:
        strains= STRAINS
    script: 'lib/collect_snps_data.py'

rule syn_snps_table:
    input:
        snps_list=temp(os.path.join(TMP_D, 'DNA_Pool1_final.tab')),
        ref_gbk=config['reference_annotation']
    output:
        syn_output= os.path.join(TMP_D, 'syn_SNPs.tab')
    params:
        script_f='lib/snps/Snp2Amino.edit.py'
    shell:
        """
        source activate py27
        python {params.script_f} -f {input.snps_list} -g {input.ref_gbk} \
        -n syn \
        -o {output.syn_output} 
        source deactivate
        """

rule nonsyn_snps_table:
    input:
        snps_list=temp(os.path.join(TMP_D, 'DNA_Pool1_final.tab')),
        ref_gbk=config['reference_annotation']
    output: 
        non_syn_output= os.path.join(TMP_D, 'non-syn_SNPs.tab'),
    params:
        script_f='lib/snps/Snp2Amino.edit.py'
    shell:
        """
        source activate py27
        python {params.script_f} -f {input.snps_list} -g {input.ref_gbk} \
        -n non-syn \
        -o {output.non_syn_output} 
        source deactivate
        """

rule all_snps_table:
    input:
        snps_list=temp(os.path.join(TMP_D, 'DNA_Pool1_final.tab')),
        ref_gbk=config['reference_annotation']
    output:
        all_snps_tab=os.path.join(TMP_D, 'all_SNPs.tsv')
    params:
        script_f='lib/snps/Snp2Amino.edit.py'
    shell:
        """
        source activate py27
        python {params.script_f} -f {input.snps_list} -g {input.ref_gbk} \
        -o {output.all_snps_tab} 
        source deactivate
        """

rule all_snps_list:
    input:
        ## the exact files should also be, so they cannot be deleted before 
        flatcount_file=expand(os.path.join(TMP_D, '{strain}', SOFTWARE['mapper'], 'dna_paired.flatcount'), strain=STRAINS),
        snp_vcf_file=expand(os.path.join(TMP_D, '{strain}', SOFTWARE['mapper'], 'dna.flt.vcf'), strain=STRAINS),
        flt_files=expand("{strain}.flatcount", 
            strain=STRAINS),
        flt_vcf_files=expand("{strain}.flt.vcf", 
            strain=STRAINS),
        dict_file= os.path.join(TMP_D, 'dict.txt'),
        anno_f='annotations_for_snps.tab'
    output: 
        snps_list=temp(os.path.join(TMP_D, 'DNA_Pool1_final.tab'))
    params:
        script_f='lib/snps/mutation_table.py'
    shell:
        """
        source activate py27
        python {params.script_f}  -f {input.dict_f} \
        -a {input.anno_f} \
        -o {output.snps_list}
        source deactivate
        """

rule for_tab_compute_annot_file:
    input:
        ref_gbk=config['reference_annotation']
    output:
        anno_f=temp('annotations_for_snps.tab')
    params:
        species= config['species']
    script:'lib/create_dict_for_snps.py'


rule create_dict_file: 
    input: 
        flt_vcf_files=expand("{strain}.flt.vcf", 
            strain=STRAINS),
        flatcount_files=expand("{strain}.flatcount", 
            strain=STRAINS)
    output:
        dict_file= temp(os.path.join(TMP_D, 'dict.txt'))
    params:
        strains=STRAINS
    shell:
        """
        echo {input.flt_vcf_files}| \
        sed 's/\.flt\.vcf\W*/\\n/g'| \
        grep '\w' > {output.dict_file}
        """


rule for_tab_copy_files:
    ## duplicate everyone to create 'STRAIN.flt.vcf' and 'STRAIN.flatcount'
    input: 
        flatcount_file=lambda wildcards: os.path.join(TMP_D, wildcards.strain, SOFTWARE['mapper'], 'dna_paired.flatcount'),
        snp_vcf_file=lambda wildcards: os.path.join(TMP_D, wildcards.strain, SOFTWARE['mapper'], 'dna.flt.vcf')
    output:
        flatcount=temp("{strain}.flatcount"),
        flt_vcf=temp("{strain}.flt.vcf")
    shell:
        """
        ln -s {input.snp_vcf_file} {output.flt_vcf}
        ln -s {input.flatcount_file} {output.flatcount}
        """

rule for_snps_convert_to_art:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/{mapper}/dna_paired.sam'
    output:
        ART=temp('{TMP_D}/{strain}/{mapper}/dna_paired.art')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params[SAM2ART]} -s 2 -p -4 {input.STAMPY_SAM} > {output} '

rule for_snps_convert_to_sin:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/{mapper}/dna_paired.sam'
    output:
        SIN=temp('{TMP_D}/{strain}/{mapper}/dna_paired.sin')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params[SAM2ART]} -s 2 -l -p -4 {input.STAMPY_SAM}  > {output} '

rule for_snps_convert_to_flatcount:
    input: 
        STAMPY_SAM='{TMP_D}/{strain}/{mapper}/dna_paired.sam'
    output:
        FLTCNT=temp('{TMP_D}/{strain}/{mapper}/dna_paired.flatcount')
    params:
        SAM2ART='lib/expr/sam2art.pl'
    shell:
        '{params[SAM2ART]} -f -s 2 -p {input.STAMPY_SAM}  > {output} '
