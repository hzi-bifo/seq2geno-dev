#!/usr/bin/perl -w
# art2genecount.pl
# AUTHOR: Andreas Dötsch, Denitsa Eckweiler
# LAST REVISED: May 2012
# 

## reads art file and writes gene count for DESeq
use strict;
use Getopt::Std;

#variables

my $usage = "\n\nusage: $0 [-t  <reference type> -c] -a <art file> -r <reference file>\n".
            "Reads a reference file of the specified format and and art file and write gene counts for expression analysis.\n\n".
	    "-a <art file>      \tart (artemis readable pileup) file containing information on expression of both strands. Should be created with \"sinister\" option.\n".
	    "-t <reference type>\ttype of reference file used. Currently only \"1\" allowed => gtf file.\n".
	    "-r <reference file>\tFile containing the reference genes. Format is specified with the -t option.\n".
	    "-b\t\tBlank mode. Only readcounts are printed in the output without gene IDs and anti-sense counts.\n".
	    "-c                 \tif this flag is set, only coding sequences (no ncRNA, tRNA, rRNA) are considered from the annotation.\n\n";
our($opt_a,$opt_t,$opt_r,$opt_c,$opt_b);
getopts('a:t:r:cb') or die $usage;

if(!defined($opt_a)){die $usage}
if(!defined($opt_t)){die $usage}
if(!defined($opt_r)){die $usage}
if(!defined($opt_b)){ $opt_b = 0 } else { $opt_b = 1 }
my $art_file	= $opt_a; 
my $type 	= $opt_t; 
my $ref_file	= $opt_r;
my $blankmode   = $opt_b;
my $cds_only;
if(defined($opt_c)){$cds_only = 1} else {$cds_only = 0}
my ($genomesize,$tmp,$i,$line,$genetype,$genestart,$geneend,$genestrand,$genetxt,$geneID,$n_genes); #skalar variables
my (@plus,@minus,@parts); #arrays
my (%genecount,%antisensecount); #hashes

open(ART,$art_file)
	  or die "Error reading file: $!\n";
open(REF,$ref_file)
	  or die "Error reading file: $!\n";

#$tmp = `wc -l $art_file`;
#@parts = split(/\s/,$tmp);
#$genomesize = $parts[0]; #genomesize is now included in the tab annotation file header
#print "$genomesize\n";

#### initalize vectors ####
#check for header in annotation
$line = <REF>;
if(substr($line,0,1) eq "@"){
	#header found, first line should contain genome size
	@parts = split(/\|/,$line);
	$genomesize = $parts[2];
	chomp($genomesize)
#	print STDERR "Genomesize found in annotation header: $genomesize\n";
}
for($i=1; $i<=$genomesize;$i++){
	$plus[$i] = 0; 
	$minus[$i] = 0;
}

#### parse ART ######
while(defined($line = <ART>)){
	#skip header
	if(substr($line,0,1) eq "#"){next}

	#read strand coverage	
	chomp($line);
	@parts = split(/\s/,$line);
	$i = $parts[0];
	$plus[$i] = $parts[1];
	$minus[$i] = $parts[2];
}

##### parse annotation #####
$n_genes = 0;
while(defined($line = <REF>)){
	chomp($line);
	#parse current gene
	if($type eq "gtf"){
		#use GTF format
		@parts = split(/\t/,$line);
#		$genetype      = ???, TODO
		$genestart = $parts[3];
		$geneend   = $parts[4];
		$genestrand= $parts[6];
		$genetxt   = $parts[8];
		@parts = split(/"/,$genetxt);
		$geneID    = $parts[1];
	}elsif($type eq "tab"){
		#skip header
		if(substr($line,0,1) eq "@"){next}

		#use tab separated format from pseudomonas.com
		@parts = split(/\t/,$line);
		$genetype   = $parts[4];
		$genestart  = $parts[5];
		$geneend    = $parts[6];

		$genestrand = $parts[7];
		if(defined($parts[8])){
			$geneID = "$parts[3],$parts[8]";
		}else{
			$geneID = $parts[3];
		}
	}else{
		print STDERR "No format description found!\n";
	}
	$n_genes++;

	#skip non-coding genes, if this option is set
	if($cds_only&&($genetype ne "CDS")){next}

	#check annotation for errors
	if(!( ($genestart =~ /^[0-9]+$/) & ($geneend =~ /^[0-9]+$/) )){
		print STDERR "!!! Annotation error for gene $geneID !!!\n";
	}

	#count reads
	for($i = $genestart; $i <= $geneend; $i++){
		if($i > $genomesize){
			print STDERR "Out of range values for gene $geneID\n";
			last;
		}
		if($genestrand eq "+"){
			$genecount{$geneID} += $plus[$i];
			$antisensecount{$geneID} += $minus[$i];
		}elsif($genestrand eq "-"){
			$genecount{$geneID} += $minus[$i];
			$antisensecount{$geneID} += $plus[$i];
		}else{
			print STDERR "Hmmm, inconsistent strand information for gene $geneID ($genestart:$geneend).\n";
		}
	}
}

foreach $geneID (sort keys %genecount){
	if($blankmode){
		print "$genecount{$geneID}\n";
	} else {
		print "$geneID\t$genecount{$geneID}\t$antisensecount{$geneID}\n";
	}
}

close ART;
close REF;