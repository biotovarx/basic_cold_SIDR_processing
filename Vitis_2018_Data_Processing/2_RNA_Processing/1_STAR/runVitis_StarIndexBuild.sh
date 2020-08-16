#!/bin/bash

script=$1
scriptstart=$(date "+%s")
now=$(date "+%Y-%m-%d_%H%M%S")
log="${script%.*}_log_$now.txt"


STAR --runThreadN 4 --runMode genomeGenerate --genomeDir star_grape_index/ --genomeFastaFiles ~/genomes/vitis/12Xv2_genome_and_annotation/12Xv2_grapevine_genome_assembly.fa --sjdbGTFfile ~/genomes/vitis/12Xv2_genome_and_annotation/vitis_vinifera_gene_annotation_on_V2_20.gtf --sjdbOverhang 74

scriptend=$(date "+%s")
runtime=$((($scriptend - $scriptstart)/60))


echo
echo "-----------------------------------------------------------------" >> $log
echo "-----------------------------------------------------------------" >> $log
echo "" >> $log
echo "Total Runtime (Minutes): " $runtime >> $log
echo "" >> $log
echo "-----------------------------------------------------------------" >> $log
echo "-----------------------------------------------------------------" >> $log
