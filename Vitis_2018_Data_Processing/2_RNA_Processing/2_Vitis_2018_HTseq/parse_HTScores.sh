#!/bin/bash



dataset_directory=$1

outdir="htseq_count_QC"

mkdir $outdir

for loopfile in $dataset_directory*
do

thisfile=${loopfile%.txt}
fileRoot=${thisfile##*/}
outputRoot=$fileRoot"_qc.txt"

tail -5 $loopfile > "$outdir/$outputRoot"



done


