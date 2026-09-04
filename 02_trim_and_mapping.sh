#!/usr/bin/bash


# Setting default parameter values
input_dir=$PWD
output_dir=$PWD
reference_genome="/home/sunchangbin/"
thread=10


# Reading in arguments
while getopts i:o:g:t:h opt
do
	case $opt in
	i)
		input_dir=$OPTARG
		;;
	o)
		output_dir=$OPTARG
		;;
	g)
		reference_genome=$OPTARG
		;;
	t)
		thread=$OPTARG
		;;
	h)
		echo "Usage:	trim_and_mapping.sh [-i input_dir] [-o output_dir] [-g reference_genome] [-t thread]"
		echo ""
		echo "Where:"
		echo "-i		Path to input directory containing FASTQ files [defaults to the working directory]"
		echo "-o		Path to output directory where to write merged FASTQ files [defaults to the working directory]"
		echo "-g		Path and base name to the reference genome indexed files that will be used for alignment ‌[defaults to the human GRCh38 reference genome]."
		echo "-t		Number of threads for parallel processing‌ [defaults to 10]"
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


if [[ ! -f "${reference_genome_path}.1.bt2" ]]
then
        echo "ERROR: Bowtie2 index not found."
        exit 2
fi

cd $output_dir
####make directoriers
mkdir fastqc
mkdir 01.cleanData
mkdir result

###
ls $input_dir/*gz | xargs fastqc -o fastqc/ -t 12
####use fastp to remove adapter contamination and low quality reads
for i in *_R1.fq.gz ; do fastp -i $i -I ${i%_*}_R2.fq.gz -o 01.cleanData/${i%%_*}_clean_R1.fq.gz  -O 01.cleanData/${i%%_*}_clean_R2.fq.gz --thread=16 --complexity_threshold=30 –n_base_limit=5 -l 50 -g -x -F 30 -f 3 ; done

###After adapter trimming, bowtie2 for mapping to reference genome
###samtools change sam files to bam files, sort and index bam files
### picard remove PCR duplcate reads
cd 01.cleanData

ls *clean_R2.fq.gz | sed 's/_clean_R2.fq.gz//' > list_clean

for i in `cat list_clean`
do
   echo $i
   bowtie2 -x $reference_genome -p $thread --local --dovetail --phred33 -1 "$i"_clean_R1.fq.gz -2 "$i"_clean_R2.fq.gz | samtools sort -@ 20 -O bam -o - > $output_dir/result/${i}.bam 
done

cd $output_dir/result
for i in *.bam ; do
    echo $i
    samtools index $i -@ $thread
    samtools flagstat $i -@ $thread > ${i}.flagstst
    java -jar picard MarkDuplicates.jar I=$i O=${i%.*}_dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true M=${i%.*}_dup.picarddup.txt
    samtools index ${i%.*}_dup.bam -@ $thread
    samtools flagstat ${i%.*}_dup.bam -@ $thread > ${i%.*}_dup.bam.flagstat
    samtools view -@ $thread -hb -q 30 ${i%.*}_dup.bam | samtools sort -@ $thread -o ${i%.*}_last.bam - 
    samtools index ${i%.*}_last.bam -@ $thread
    samtools flagstat ${i%.*}_last.bam -@ $thread > ${i%.*}_last.bam.flagstat
    rm ${i}.bam
done
###Remove temp files
rm *-1.bam *-2.bam *dup.bam ${i%.*}_dup.bam.bai

echo "trim_and_mapping:   ...done!"

