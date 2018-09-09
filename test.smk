
rule all:
    input:
        'pyVersion2'
 
    
rule b:
    input: 
        'test.txt'
    output: 
        'pyVersion2'
    shell:
        """
        python -V 1&2> {output}
        """
rule a:
    input:
        'test.smk'
    output:
        'test.txt',
        'pyVersion1'
    shell:
        """
        
        echo 'rule 1' > {output[0]} 
        python -V 2> {output[1]}
        """

