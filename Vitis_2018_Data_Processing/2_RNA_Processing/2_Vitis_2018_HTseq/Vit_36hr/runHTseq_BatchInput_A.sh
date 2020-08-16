#!/bin/bash

# This script will run the HTseq-count on a series of bam files
# The input is a directory of STAR output directories

# MAKE SURE YOU HAVE HTseq IN THE PATH; CHANGE CONDA ENV IF NECESSARY PRIOR TO SCRIPT!

scriptstart=$(date "+%s")
# Clear the terminal screen and move down one line:
clear
echo ""


# ***************************************************************************
# › › › › CHECK ARGUMENTS
# ***************************************************************************
# Check for the required list of files in the arguments:
if [ $# -lt 1 ]
then
echo "**********************************************"
echo "*** ERROR:  NO FILES SPECIFIED IN ARGUMENT ***"
echo "Please specify the directory with the sequence files as the argument"
echo "when calling the script \"$0\""
echo "Example:  \"sh $0 *.fastq\""
echo "**********************************************"
echo
exit 1
fi


# Display the number of files and their names:
echo "Number of files passed as a set in one argument = $#"
echo "_____________________________________________________"
echo
echo "Files:"
echo
for i in $@
do
echo $i
done
echo "_____________________________________________________"
echo

# ***************************************************************************
#›› › › SET VARIABLES
# ***************************************************************************
now=$(date "+%Y-%m-%d_%H%M%S")
script=$0


# Constants
HTSEQ_OUTPUT="~/data_processing/htseq_counting/Vit_36hr/Vitis_STAR_HTseq_Output_$now"
GENE_FILE="~/genomes/vitis/12Xv2_genome_and_annotation/vitis_vinifera_gene_annotation_on_V2_20.gtf"

# ***************************************************************************
#›› › › SET UP LOG
#›› › › FILES
# ***************************************************************************

mkdir $HTSEQ_OUTPUT

log="$HTSEQ_OUTPUT/${script%.*}_log_$now.txt"


# Header for log:
echo "_____________________________________________________" > $log
echo >> $log
echo "Log of STAR run `date`:" >> $log
echo "_____________________________________________________" >> $log
echo >> $log
echo "Number of files passed as a set in one argument = $#" >> $log
echo "_____________________________________________________" >> $log
echo >> $log
echo "Files:" >> $log
echo >> $log
for i in $@
do
echo $i >> $log
done
echo "_____________________________________________________" >> $log
echo >> $log
echo "\$0 = $0" >> $log
echo "\$now = $now" >> $log
echo "\$script = $script" >> $log
echo "\$log = $log" >> $log
echo "_____________________________________________________" >> $log
echo >> $log


# ***************************************************************************
#				RUN HTseq-count
#
# ***************************************************************************

for loopDir in $@
do

	if [ ! -d $loopDir ]; then
		continue
	fi


	# Unzip file using gunzip and adjust filename accordingly:
	#gunzip $loopfile

	dirroot=${loopDir##*/}
	fileroot=${dirroot%_STR*}
	infile=$loopDir"/"$fileroot".bam"
	outRoot=${fileroot%_TR*}
	thisOutCount=$outRoot"_counts.txt"

	echo "dirroot = $dirroot" >> $log
	echo "fileroot = $fileroot" >> $log
	echo "infile = $infile" >> $log
	echo "thisOutCount = $thisOutCount" >> $log


	htseq-count -f bam --stranded=no $infile $GENE_FILE > $HTSEQ_OUTPUT"/"$thisOutCount




done


