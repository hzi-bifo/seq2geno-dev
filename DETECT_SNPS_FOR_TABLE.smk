'''
Purpose:
Detect variants for the snps table.

Output:
        vcf=temp("{TMP_D}/{strain}/{mapper}/dna.flt.vcf")
'''
    
rule for_tab_filter_vcf:
    input:
        bcf="{TMP_D}/{strain}/{mapper}/dna.raw.bcf"
    output:
        vcf=temp("{TMP_D}/{strain}/{mapper}/dna.flt.vcf")
    params: 
        VCFUTIL_EXE='lib/vcfutils.pl',
        minDepth= 0
    shell:
        """
        source activate Ariane_dna
        bcftools view {input.bcf_out} |\
        {params.VCFUTIL_EXE} varFilter -d {params.minDepth} > {output.vcf_out}
        """ 
    
rule for_tab_create_vcf:
    input:
        REF=REF_FA,
        REF_FA_INDEX=REF_FA+".fai",
        BAM="{TMP_D}/{strain}/{mapper}/dna_paired_sorted.bam",
        BAM_INDEX="{TMP_D}/{strain}/{mapper}/dna_paired_sorted.bam.bai"
    output:
        bcf=temp("{TMP_D}/{strain}/{mapper}/dna.raw.bcf")
    params: 
        minDepth= 0,
        CORES=CORES
    shell:
        """
        source activate Ariane_dna
        samtools mpileup -uf {input.REF} {input.BAM} |\
        bcftools view -bvcg - > {output.bcf_out}
        source deactivate
        """ 

rule for_tab_sort_bam:
    input:
        paired_bam="{TMP_D}/{strain}/{mapper}/dna_paired.bam"
    output:
        sorted_bam=temp("{TMP_D}/{strain}/{mapper}/dna_paired_sorted.bam"),
        sorted_bam_index=temp("{TMP_D}/{strain}/{mapper}/dna_paired_sorted.bam.bai")
    params:
        output_prefix= lambda wildcards: os.path.join(wildcards.TMP_D, wildcards.strain, wildcards.mapper, 'dna_paired_sorted')
    shell:
        """
        source activate Ariane_dna
        samtools sort {input.paired_bam} {params.output_prefix}
        samtools index {output.sorted_bam}
        source deactivate
        """


rule for_tab_sam2bam:
    input:
        STAMPY_SAM= temp('{TMP_D}/{strain}/{mapper}/dna_paired.sam')
    output:
        STAMPY_BAM= temp('{TMP_D}/{strain}/{mapper}/dna_paired.bam')
    params:
        CORES=CORES
    shell:
        """
        source activate Ariane_dna
        samtools view -bS -@ {params.CORES} \
        {input.STAMPY_SAM} > {output.STAMPY_BAM}
        source deactivate
        """

rule for_tab_stampy_mapping:
    input:
        FQ1=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads1'],
        FQ2=lambda wildcards: SAMPLES_DF.loc[wildcards.strain, 'reads2'],
        REF=REF_FA,
        REF_BWA_INDEX=REF_FA+".fai",
        REF_BWA_INDEXDICT=REF_FA+".bwt"
    output:
        STAMPY_SAM= temp('{TMP_D}/{strain}/{mapper}/dna_paired.sam')
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY_BIN=STAMPY_EXE
    shell:
        """
        source activate Ariane_dna
        {params.STAMPY_BIN} --bwaoptions=\"-q10 {input.REF}\" \
        -g {params.REF_PREFIX} -h {params.REF_PREFIX} -t{params.CORES} \
        -M {input.FQ1} {input.FQ2} > {output[STAMPY_SAM]} 
        source deactivate
        """
