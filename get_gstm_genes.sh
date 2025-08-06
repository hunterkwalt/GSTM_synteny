#!/bin/bash

#this script retrieves flanking genes from a list of genes of interest and makes a tsv file
## The tsv file contains information for making synteny plots through gggenomes

### This should be run in a directory containing folders by: 
###	1.) the name of your taxa of interest 
###	2.) a single annotation file within each directory
###	3.) a tab-separated file with your taxa of interest in one column, and the genes of interest in the second.

#working directory
wd=$PWD

#file with taxa and gene names
toi=$wd/rodent_gstm_genes.txt

#number of flanking genes on each side
bf=2
af=2

#output directory
op=$wd/synteny_output
mkdir -p $op #make output directory

#loop to pull genes
while read -r line
do
	#grep $g $file | cut -f 1
	spp=$(echo $line | cut -f 1 -d " ")
	goi=$(echo $line | cut -f 2 -d " ")
	gff=$(echo $wd/${spp}/*)
	#echo $common
	awk '$3=="gene"' $gff | grep "biotype=protein_coding" | grep -w ${goi} -B $bf -A $af | cut -f 1,4,5,7,9 
done<$toi > $op/gene_flank_locs.tmp

### I had to do Phyllotis and neotoma a little differently
po=phyllotis_neotoma.tsv
grep Phyllotis $toi > $po #subset phyllotis
grep Neotoma $toi >> $po

#loop for phyllotis
while read -r line
do
	#grep $g $file | cut -f 1
	spp=$(echo $line | cut -f 1 -d " ")
	goi=$(echo $line | cut -f 2 -d " ")
	gff=$(echo $wd/${spp}/*)
	#echo $common
	awk '$3=="gene"' $gff | grep -w ${goi} -B $bf -A $af | cut -f 1,4,5,7,9
done<$po >> $op/gene_flank_locs.tmp

#get rid of redundant genes
sort -u $op/gene_flank_locs.tmp > $op/gene_flank_locs.tsv

cut -f 5 $op/gene_flank_locs.tsv | cut -f 1 -d ";" | sed "s/ID=//" | sed "s/gene-//" | sed "s/gene_id//" | sed 's/"//g' | sed "s/ //g" > $op/names.tmp

paste $op/gene_flank_locs.tsv $op/names.tmp > gene_flank_locs_final.tsv

rm $op/gene_flank_locs.tmp #delete file with redundant genes
rm $op/names.tmp # delete names



 
 
 
