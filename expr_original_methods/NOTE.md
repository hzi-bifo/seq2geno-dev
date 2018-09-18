my_stampy_pipeline CH2522 CH2522.fastq /data3/reference_sequences/stampy-references/new_stampy-references/Pseudomonas_aeruginosa_PA14.fasta /data3/reference_sequences/Pseudomonas_aeruginosa_PA14_annotation_with_ncRNAs_07_2011_12genes.tab /data3/reference_sequences/Pseudomonas_aeruginosa_PA14_12genes_R_annotation 2> CH2522.log&
art2gene_coverage.pl -t tab -a CH2522.art -r /data3/reference_sequences/Pseudomonas_aeruginosa_PA14_annotation_with_ncRNAs_07_2011.tab  >  CH2522_coverage&

	my_stampy_pipeline
	    - sam2art.pl 
		-s 2 $prefix.sam > $prefix.art
		-s 2 -l $prefix.sam > $prefix.sin
		-f -s 2 $prefix.sam > $prefix.flatcount
	    - art2genecount.pl 
		-b -a $prefix.sin -t tab -r $annotation > $prefix.rpg
		- annotation: path to annotation file for R
	    - sam_statistics.pl 
		-r $prefix.sam > $prefix.stats
	    - genes_statistics.R 
		$prefix.rpg $Rannotation $prefix.stats $prefix
	art2cov.py 
		-f dictionary -s SE -c 1 -o coverage_cut1.txt
		- dictionary: prefix of ".art" and ".rstats"
