#!/bin/bash

clear
echo

##### Calulate runtime #########
scriptstart=$(date "+%s")

echo "This script will run FastQC on the requested files."
echo
echo

# ***************************************************************************
#	SET VARIABLES
# ***************************************************************************
now=$(date "+%Y-%m-%d_%H%M%S")
script=$0
dataset_directory=$1
species="Athaliana"

logfile="$PWD/${script%.*}_$species_log_$now.txt"
callQC="/usr/local/bin/FastQC_11.7/fastqc"
dirout="$PWD/${script%.*}_$species_output_$now"

#*************************************************************
#	 SET UP LOG
#*************************************************************

echo >> $logfile
echo "Log of FastQC run `date`:" >> $logfile
echo >> $logfile
echo >> $logfile


# ***************************************************************************
#	 RUN FASTQC
# ***************************************************************************
mkdir $dirout


for file in $dataset_directory*
do

echo "FastX-CL processing file $file - `date`" >> $logfile

$callQC -a ~/data_processing/Preprocessing/TruSeq_Adapter.txt $file --outdir=$dirout

echo "--------------------------------------------">> $logfile
echo >> $logfile
done

scriptend=$(date "+%s")
runtime=$((($scriptend - $scriptstart)/60))


echo
echo "-----------------------------------------------------------------" >> $logfile
echo "-----------------------------------------------------------------" >> $logfile
echo "" >> $logfile
echo "Total Runtime (Minutes): " $runtime >> $logfile
echo "" >> $logfile
echo "-----------------------------------------------------------------" >> $logfile
echo "-----------------------------------------------------------------" >> $logfile


