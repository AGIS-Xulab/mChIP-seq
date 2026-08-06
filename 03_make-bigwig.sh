#!/usr/bin/bash

# Setting default parameter values
input_dir=$PWD
output_dir=$PWD
effectiveGenomeSize=3042914375
thread=10


# Reading in arguments
while getopts i:o:s:t:h opt
do
	case $opt in
	i)
		input_dir=$OPTARG
		;;
	o)
		output_dir=$OPTARG
		;;
	s)
		effectiveGenomeSize=$OPTARG
		;;
	t)
		thread=$OPTARG
		;;
	h)
		echo "Usage:	make-bigwig.sh [-i input_dir] [-o output_dir] [-s effectiveGenomeSize] [-t thread]"
		echo ""
		echo "Where:"
		echo "-i		Path to input directory containing bedGraph files with methylation calls [defaults to the working directory]"
		echo "-o		Path to output directory where to store compressed bigwig files for visualisation [defaults to the working directory]"
		echo "-s		Effecitive genome size. [defaults to size of GRCh38]"
		echo "-t		Number of threads for parallel processing‌ [defaults to 10]"
		echo ""
		exit 1
		;;
	esac
done

# Validating arguments
echo "[mapping-stats]:	Validating arguments..."

if [[ ! -d $input_dir ]]
then
		echo "[mapping-stats]:	ERROR: Input directory not found."
		exit 2
fi 

if [[ ! -d $output_dir ]]
then
		echo "[mapping-stats]:	ERROR: Output directory not found."
        exit 2
fi 


for i in *_dup.bam ;do

         bamCoverage -v -p $thread -b $i -o ${i%.*}.bw -bs 10 --effectiveGenomeSize $effectiveGenomeSize --normalizeUsing RPGC

done

echo "[make-bigwig]:       ...done!"
