#!/bin/bash

script=$0
scriptstart=$(date "+%s")
now=$(date "+%Y-%m-%d_%H%M%S")
log="${script%.*}_log_$now.txt"


STAR --runThreadN 4 --runMode genomeGenerate --genomeDir star_AT_index --genomeFastaFiles ~/genomes/AT/Tair10_Genome/AT_Tair10_genome.fa --sjdbGTFfile ~/genomes/AT/Tair10_Genome/Araport11_GFF3_genes_transposons.201606.gtf --sjdbOverhang 73

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
