#!/bin/bash

# This script will run the STAR aligner on a series of fastq files


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
for i in $1*
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
dataset_directory=$1
datasetName="AF_TR_QT_QF_QF"

# Constants
STARINDEX="~/data_processing/star_alignment/Arabidopsis_2018_STAR/star_AT_index"
OPTIONS="--runThreadN 6"


# Used for output to Cufflinks; NOT needed for 3'prime sequence
#CUFFOPTS="--outSAMstrandField intronMotif"

# These options were used by Karl Kremling for his maize paper, in the Buckler lab
THREEprimeOPS="--outFilterMultimapNmax 10 --outFilterMismatchNoverLmax 0.04 --outFilterIntronMotifs RemoveNoncanonicalUnannotated"
STAR_OUTPUT="Arabidopsis_"$datasetName"_2018_STAR_Genome_$now"


# ***************************************************************************
#›› › › SET UP LOG
#›› › › FILES
# ***************************************************************************

mkdir $STAR_OUTPUT
cd $STAR_OUTPUT

log="$PWD/${script%.*}_log_$now.txt"
alignSumLog="$PWD/${script%.*}_alignment_summary_$now.txt"

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
for i in $dataset_directory*
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

# Header for alignment summary file:

echo "_____________________________________________________" > $alignSumLog
echo >> $alignSumLog
echo "Log of STAR run `date`:" >> $alignSumLog
echo "_____________________________________________________" >> $alignSumLog
echo >> $alignSumLog
echo "Number of files passed as a set in one argument = $#" >> $alignSumLog
echo "_____________________________________________________" >> $alignSumLog
echo >> $alignSumLog
echo "Files:" >> $alignSumLog
echo >> $alignSumLog
for i in $dataset_directory*
do
echo $i >> $alignSumLog
done
echo "_____________________________________________________" >> $alignSumLog
echo >> $alignSumLog
echo "\$0 = $0" >> $alignSumLog
echo "\$now = $now" >> $alignSumLog
echo "\$script = $script" >> $alignSumLog
echo "\$log = $log" >> $alignSumLog
echo "_____________________________________________________" >> $alignSumLog
echo >> $alignSumLog


# ***************************************************************************
#›› › › RUN STAR
#›› › › ALIGNMENT
# ***************************************************************************

for loopfile in $dataset_directory*
do
# Unzip file using gunzip and adjust filename accordingly:
#gunzip $loopfile
thisfile=${loopfile%.fastq.gz}
fileRoot=${thisfile##*/}
outputRoot=${fileRoot%_CL*}
thisOutDir=$outputRoot"_STRout"
starout=$outputRoot"_"
mkdir $PWD"/"$thisOutDir
cd $thisOutDir

thisCall="STAR $OPTIONS $THREEprimeOPS --outReadsUnmapped Fastx --genomeDir $STARINDEX --readFilesIn $loopfile --readFilesCommand gunzip -c --outFileNamePrefix $starout"

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Record date, time, and command call to log file:
echo "_________________________________________________________________" >> $log
echo >> $log
echo `date` >> $log
echo $thisCall >> $log
echo "This file: " $thisfile >> $log
echo "File root: " $fileRoot >> $log
echo "Output Root " $outputRoot >> $log
echo "This out dir: " $thisOutDir >> $log
echo "-----------------------------------------------------------------" >> $log

# Execute command call:
$thisCall

# Add alignment statistics to alignment summary file:
summaryFile=$starout"Log.final.out"
echo "_________________________________________________________________" >> $alignSumLog
echo >> $alignSumLog
echo `date` >> $alignSumLog
echo "Output Root " $outputRoot >> $alignSumLog
echo "This out dir: " $thisOutDir >> $alignSumLog
echo >> $alignSumLog
cat $summaryFile >> $alignSumLog
echo >> $alignSumLog
echo "-----------------------------------------------------------------" >> $alignSumLog


# Convert the SAM output file to a BAM file

samOptions="sort -n -@ 6"

samtools $samOptions -o $outputRoot".bam" $outputRoot"_Aligned.out.sam"


gzip $outputRoot"_Aligned.out.sam"



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++›
cd ".."
# Zip original input file ($thisfile) and final output file ($outfile)
# using gzip:
#gzip $thisfile
done

scriptend=$(date "+%s")
runtime=$((($scriptend - $scriptstart)/60))

echo ""
echo "-----------------------------------------------------------------" >> $log
echo "-----------------------------------------------------------------" >> $log
echo "" >> $log
echo "Total Runtime (Minutes): " $runtime >> $log
echo "" >> $log
echo "-----------------------------------------------------------------" >> $log
echo "-----------------------------------------------------------------" >> $log
