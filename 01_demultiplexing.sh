#!/usr/bin/bash


# Setting default parameter values
###raw fastq files (R1 and R2) in input_dir, PE sequencing 
input_dir=$PWD
###store fastq files after demultiplexing
output_dir=$PWD

# Reading in arguments
while getopts i:o:s:h opt
do
	case $opt in
	i)
		input_dir=$OPTARG
		;;
	o)
		output_dir=$OPTARG
		;;
	s)
		barcode_file=$OPTARG
		;;
	h)
		echo "Usage:	demultiplexing.sh [-i input_dir] [-o output_dir] [-b barcode_file]"
		echo ""
		echo "Where:"
		echo "-i		Path to input directory containing FASTQ files [defaults to the working directory]"
		echo "-o		Path to output directory where to write merged FASTQ files [defaults to the working directory]"
		echo "-b		Path to a barcode file (one sample per line), ‌formatted according to fastq-multx specifications‌."
		echo ""
		exit 1
		;;
	esac
done


# Validating arguments
echo "[demultiplexing.sh]:	Validating arguments..."

if [[ ! -d $input_dir ]]
then
        echo "ERROR: Input directory not found."
        exit 2
fi 

if [[ ! -d $output_dir ]]
then
        echo "ERROR: Output directory not found."
        exit 2
fi 


if [[ ! -f $barcode_file ]]
then
        echo "ERROR: Barcode file not found"
        exit 2
fi

###demultiplexing fastq files
###Sequences of BC in Read2 5'end were used
###barcode file includes BC sequences that should be complementary reverse to sequences of Adapter
echo "Processing demultiplexing..."
fastq-multx -B ${barcode_file} ${input_dir}/*_R2.fq.gz ${input_dir}/*_R1.fq.gz -m 0 -b -x -o ${output_dir}/%_R2.fastq -o ${output_dir}/%_R1.fastq > ${output_dir}/fastq-multx.info 2>&1

###Compressing fastq files
echo "Processing fastq files..."
gzip ${output_dir}/*fastq


echo "demultiplexing:   ...done!"

