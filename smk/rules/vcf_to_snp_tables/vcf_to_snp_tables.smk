rule vcf_to_snp_tables:
    input:
        coding_vcf_gz=os.path.join(TMP_D, 'freebayes',
'multisample.vcf.coding.gz'),
        igr_vcf_gz= os.path.join(TMP_D, 'freebayes',
'multisample.vcf.igr.gz'),
        ref_gbk=config['reference_annotation']
    output:
        syn_snps_table= config['syn_snps_table'],
        nonsyn_snps_table= config['nonsyn_snps_table']
    params:
        strains= DNA_READS.index.values.tolist(),
        bcftools_bin= 'bcftools'
    script: 'vcf_to_feat_table.py'
