'''
Purpose:
Detect variants for the snps table.

Output:
        vcf=temp("{TMP_D,"^[^/]+$"}/{strain}/samtools/dna.flt.vcf")
'''
    
rule for_tab_filter_vcf:
    input:
        bcf='{TMP_D}/{strain}/samtools/tab_raw.bcf'
    output:
        vcf='{TMP_D}/{strain}/samtools/dna.flt.vcf'
    params: 
        VCFUTIL_EXE='lib/vcfutils.pl',
        minDepth= 0
    shell:
        """
        source activate Ariane_dna
        bcftools view {input.bcf} |\
        {params.VCFUTIL_EXE} varFilter -d {params.minDepth} > {output.vcf}
        """ 
    
rule for_tab_create_vcf:
    input:
        REF=REF_FA,
        REF_FA_INDEX=REF_FA+".fai",
        BAM='{TMP_D}/{strain}/stampy/tab_dna_sorted.bam',
        BAM_INDEX='{TMP_D}/{strain}/stampy/tab_dna_sorted.bam.bai'
    output:
        bcf=temp('{TMP_D}/{strain}/samtools/tab_raw.bcf')
    params: 
        minDepth= 0,
        CORES=CORES
    shell:
        """
        source activate Ariane_dna
        samtools mpileup -uf {input.REF} {input.BAM} |\
        bcftools view -bvcg - > {output.bcf}
        source deactivate
        """ 

rule for_tab_sort_bam:
    input:
        paired_bam='{TMP_D}/{strain}/stampy/tab_dna.bam'
    output:
        sorted_bam=temp("{TMP_D}/{strain}/stampy/tab_dna_sorted.bam"),
        sorted_bam_index=temp("{TMP_D}/{strain}/stampy/tab_dna_sorted.bam.bai")
    params:
        output_prefix= lambda wildcards: os.path.join(wildcards.TMP_D,
wildcards.strain, 'stampy', 'tab_dna_sorted')
    shell:
        """
        source activate Ariane_dna
        samtools sort {input.paired_bam} {params.output_prefix}
        samtools index {output.sorted_bam}
        source deactivate
        """


rule for_tab_sam2bam:
    input:
        STAMPY_SAM= temp('{TMP_D}/{strain}/stampy/tab_dna.sam')
    output:
        STAMPY_BAM= temp('{TMP_D}/{strain}/stampy/tab_dna.bam')
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
        STAMPY_SAM= temp('{TMP_D}/{strain}/stampy/tab_dna.sam')
    params:
        CORES=CORES,
        REF_PREFIX=REF_FA,
        STAMPY_BIN=STAMPY_EXE
    shell:
        """
        source activate Ariane_dna
        {params.STAMPY_BIN} --bwaoptions=\"-q10 {input.REF}\" \
-g {params.REF_PREFIX} \
-h {params.REF_PREFIX} \
-M {input.FQ1} {input.FQ2} > {output.STAMPY_SAM} 
        source deactivate
        """
