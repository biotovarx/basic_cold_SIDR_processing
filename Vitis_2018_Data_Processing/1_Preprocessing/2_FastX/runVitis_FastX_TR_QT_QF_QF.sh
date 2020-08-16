#!/bin/bash

# Run Fastx processing on a series of fastq files:


##### Calulate runtime #########
scriptstart=$(date "+%s")


# Clear the terminal screen and move down one line:
clear
echo


# *******************************************************
#		Set Variables
# ******************************************************

# Set date
now=$(date "+%Y-%m-%d_%H%M%S")

# Set program paths
callQC="/usr/local/bin/FastQC_11.7/fastqc"
callTR="/usr/local/bin/FastX/0.0.13/fastx_trimmer"
callQT="/usr/local/bin/FastX/0.0.13/fastq_quality_trimmer"
callQF="/usr/local/bin/FastX/0.0.13/fastq_quality_filter"

# Set FastX options
options="-v -Q33"

# Change based on dataset
species="Vitis"

# Set script variables
script=$0
dataset_directory=$1

scriptDir=${script%/*}
echo "$scriptDir"


 # Use these options when calling the script from another directory
fastqcDir="${script%.*}_"$species"_fastQC_output_$now"
mkdir $fastqcDir
log="${script%.*}_"$species"_log_$now.txt"
fastxOutDir="${script%.*}_"$species"_fastX_output_$now"
mkdir $fastxOutDir


countsfile="$fastqcDir/fastq_read_counts_$now.csv"
lengthsfile="$fastqcDir/fastq_read_lengths_$now.csv"


# ***************************************************************************
# 				CHECK ARGUMENTS
# ***************************************************************************
# Check for the required list of files in the arguments:
if [ $# -lt 1 ]
then
	echo "**********************************************"
	echo "*** ERROR:  NO FILES SPECIFIED IN ARGUMENT ***"
	echo "Please specify a list of files as the argument"
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
for i in $dataset_directory*
do
	echo $i
done
echo "_____________________________________________________"
echo

# ***************************************************************************
#				SET UP LOG
#				FILES
# ***************************************************************************
# Header for log:
echo "__________________________________________________" > $log
echo >> $log
echo "Log of FastX Quality Trimmer run `date`:" >> $log
echo "__________________________________________________" >> $log
echo >> $log

# Output the number of files and their names to log file:
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


# Header for countsfile:
echo "Filename,Process,countRead" > $countsfile

# Header for lengthsfile:
# Set the length of Array equal to the max length of reads
lengthsArr=("Filename" {1..90})

# "IFS" equals internal field separator; used to create ',' separated lines
lengthsStr=$( IFS=$','; echo "${lengthsArr[*]}" )
echo $lengthsStr > $lengthsfile

# ***************************************************************************
#				RUN FASTX
#				PROCESSES
# ***************************************************************************
for loopfile in $dataset_directory*
do
	# Unzip file using gunzip and adjust filename accordingly:
	gunzip $loopfile
	thisfile=${loopfile%.gz}


	## Use these options when not unzipping file
        #thisfile=$loopfile

	# Count reads in the initial input file:
	read linecount filename <<< $(wc -l $thisfile)
	echo "${filename##*/},Initial,$((linecount / 4))" >> $countsfile


	# Set initial value for file name:
	outfile=$thisfile

	# Trim based from either 5' or 3' end, based on optFX parameters:
	# ------------------------------------------------------------
	# FastX settings:
	#	"TR" / $callTR -> Trimmer

	process="TR"
	callFX=$callTR
	optFX="-f 13"
	# ------------------------------------------------------------
	infile=$outfile
	outfile=${infile%.*}"_$process.fastq"
	echo "FastX-$process processing file $infile - `date`" >> $log
	$callFX $options $optFX -i $infile -o $outfile >> $log
	# ------------------------------------------------------------
	echo "__________________________________________________" >> $log
	echo >> $log

	# Count reads in this output file:
	read linecount filename <<< $(wc -l $outfile)
        echo "${filename##*/},Post-$process,$((linecount / 4))" >> $countsfile

	# Trim low quality sequence from 3' ends of reads:
	# ------------------------------------------------------------
	# FastX settings:
	#	"QT" / $callQT -> Quality Trimmer

	process="QT"
	callFX=$callQT
	optFX="-t 20 -l 38"
	# ------------------------------------------------------------
	infile=$outfile
	outfile=${infile%.*}"_$process.fastq"
	echo "FastX-$process processing file $infile - `date`" >> $log
	$callFX $options $optFX -i $infile -o $outfile >> $log
	# ------------------------------------------------------------
	echo "__________________________________________________" >> $log
	echo >> $log
	
	# Count reads in this output file:
	read linecount filename <<< $(wc -l $outfile)
        echo "${filename##*/},Post-$process,$((linecount / 4))" >> $countsfile

	# Remove unneeded intermediate files:
	rm $infile

	# Quality filtering, part 1:
	# ------------------------------------------------------------
	# FastX settings:
	#	"QF" / $callQF -> Quality Filter

	process="QF"
	callFX=$callQF
	optFX="-q 20 -p 80"
	# ------------------------------------------------------------
	infile=$outfile
	outfile=${infile%.*}"_$process.fastq"
	echo "FastX-$process processing file $infile - `date`" >> $log
	$callFX $options $optFX -i $infile -o $outfile >> $log
	# ------------------------------------------------------------
	echo "__________________________________________________" >> $log
	echo >> $log

	# Count reads in this output file:
	read linecount filename <<< $(wc -l $outfile)
        echo "${filename##*/},Post-$process,$((linecount / 4))" >> $countsfile


	# Remove unneeded intermediate files:
	rm $infile

	# Quality filtering, part 2:
	# ------------------------------------------------------------
	# FastX settings:
	#	"QF" / $callQF -> Quality Filter

	process="QF"
	callFX=$callQF
	optFX="-q 13 -p 100"
	# ------------------------------------------------------------
	infile=$outfile
	outfile=${infile%.*}"_$process.fastq"
	echo "FastX-$process processing file $infile - `date`" >> $log
	$callFX $options $optFX -i $infile -o $outfile >> $log
	# ------------------------------------------------------------
	echo "__________________________________________________" >> $log
	echo >> $log

	# Count reads in this output file:
	read linecount filename <<< $(wc -l $outfile)
        echo "${filename##*/},Post-$process,$((linecount / 4))" >> $countsfile

	# Remove unneeded intermediate files:
	rm $infile

# ------------------------------------------------------------
	# Compile read lengths distribution for this final output file and save
	# results as a row in the csv file ($lengthsfile):
	for i in  ${!lengthsArr[*]}
	do
		lengthsArr[$i]=0
	done
	lengthsArr[0]=$outfile
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	i=0
	j=0
	firstID=0
	# Reads one line at a time;
	# -r -> backslash does not act as an escape character.
	# The backslash is considered to be part of the line.
	# In particular, a backslash-newline pair may not be used as a line continuation

	while read -r line
	do
		((j++))
		# Checks if firstID is 0 and if line starts with a '@'
		if [[ $firstID == 0 && $line == @* ]]
		then
			firstID=1
			# resets j to 3 for so that the next line (read line) j will equal 4
			j=3
		fi
		if (( $j % 4 == 0))
		then
			# Incriments the length array at the index equal to the length of the line (length of read)
			((lengthsArr[${#line}]++))
		fi
	done < "$outfile"
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	lengthsStr=$( IFS=$','; echo "${lengthsArr[*]}" )
	echo $lengthsStr >> $lengthsfile
	# ------------------------------------------------------------


	# Run FastQC on the final output file:
	$callQC $outfile --outdir=$fastqcDir

	# Zip original input file ($thisfile) and final output file ($outfile)
	# using gzip:
	gzip $thisfile $outfile

	# Move fastx processing output file to fastx output directory:
	mv $outfile".gz" $fastxOutDir
done


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





