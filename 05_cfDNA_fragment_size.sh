#!/usr/bin/bash

#######################################################################
########################## filter bam files ###########################

# Define project directory
PROJECT_DIR= "path to directory"
##PROJECT_DIR=/user/project
cd "$PROJECT_DIR"

# Define/make input/output directories
bam_dir="${PROJECT_DIR}/bam_files" #Aligned BAM files
fbam_dir="${PROJECT_DIR}/01-filter_bam" #Filtered BAM files
stat_dir="${PROJECT_DIR}/stat_reports" #Alignment statistics
mkdir -p "$bam_dir" "$fbam_dir" "$stat_dir" #Create directories if they do not exist

# Obtain file path for all BAM files
bam_files=$(find $bam_dir -maxdepth 1 -name '*.bam')

# Iterate through each BAM file in directory
for bam_file in $bam_files; do
  fname=$(basename $bam_file)
  out_file="${fbam_dir}/${fname%.bam}_f1.bam"
  echo "Output file will be written to: $out_file"

  # Use samtools view to apply filtering criteria and remove chrom Y and MT
  chrom="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"
  samtools view -@ 100 -bh -f 2 -F 3844 -q 30 $bam_file ${chrom}> $out_file
  # Summary statistics for the filtered BAM file
  stat_file="${stat_dir}/${fname%.bam}_flagstat.txt"
  echo "Generating stat report: $stat_file"
  samtools flagstat $out_file > $stat_file

done
echo " Filtering bam files......Finished "


#######################################################################
########################### bam to bed ################################

# Define directories for output files
bedpe_dir="${PROJECT_DIR}/02-bamtobed" #BEDPE files
mkdir -p "$bedpe_dir" #Create directories if they do not exist


# Convert BAM files BEDPE files using bedtools
for f in $(find $fbam_dir -maxdepth 1 -iname "*.bam" -type f)
do
  id=$(basename -a -s f1.bam $f | cut -d "_" -f 1)
  echo "id:" $id
  # Files must be sorted bam before converting to BEDPE
  samtools sort -n -@ 100 $f | \
  # Convert to BEDPE file
  bedtools bamtobed -i - -bedpe |\
  # Calculate insert size and add 'chr' to chromosome field
  awk 'OFS = "\t" {print $1=$1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $6-$2}' - > $bedpe_dir/${id}.bedpe
done
echo "Converting BAM files to BEDPE files ......Finished "


#######################################################################
########################### GC content ################################

# Define directories for output files
gc_dir="${PROJECT_DIR}/03-filter_frags_gc" #filtered fragments
# path to genome
genome="/disk1/sunchangbin/genome/GRCh38/genecode/GRCh38.primary_assembly.genome.fa"
#path to blacklisted bed file
filter="/disk1/sunchangbin/prj/hg38-blacklist.v2.chr.bed"
#path to genomic gaps bed file
gap="/disk1/sunchangbin/prj/mChIP-cf-data_analysis/fragment/gaps.hg38.bed"
#Create directories if needed
mkdir -p "$out_dir" 

# Iterate through all BED files
for f in $(find $bedpe_dir -maxdepth 1 -iname "*.bedpe" -type f)
do
  # Get sample ID
  id=$(basename -a -s .bedpe $f)
  echo "id:" $id
  # BED file format: fragment start/end interval/QNAME/frag size
  awk 'OFS = "\t" {print $1, $2, $6, $7, $11}' $f |\
  # Sort by chromosome and start position
  sort -k1,1 -k2,2n |\
  # Filter blacklisted/gap regions
  bedtools subtract -a - -b $filter -A |\
  bedtools subtract -a - -b $gap -A |\
  # Obtain fragment GC content
  bedtools nuc -fi $genome -bed - | \
  # Add GC content($7)
  awk 'OFS = "\t" {print $1, $2, $3, $4, $5, $7}' - > $gc_dir/${id}_frags_gc.bed
  echo $f Complete
done
echo " Filtering and calculating fragment GC content.......Finished "


#######################################################################
######### unique and de-duplicated fragments in 5 MB bins #############

# Define directories foroutput files
frags_dir="${PROJECT_DIR}/04-bins5mb"
bins5mb="/disk1/sunchangbin/prj/mChIP-cf-data_analysis/fragment/hg38.bins5mb_filtered.bed"
mkdir -p "frags_dir" #Create directories if they do not exist


#  Iterate through all bed files
for f in $(find $gc_dir -maxdepth 1 -iname "*.bed" -type f)
do
  # Obtain ID from filename
  id=$(basename -a -s _frags_gc.bed $f)
  echo $id
  echo $f

  # Intersect fragments with 5Mb bins and unique intersections
  bedtools intersect -a $bins5mb -b $f -wo |\
  # Moved base overlap ($11) to $5 to avoid issues with uniq syntax
  awk 'OFS="\t" {ov=$11; print $1, $2, $3, $4, ov, $5, $6, $7, $8, $9, $10}' - |\
  # Sort by QNAME ($9) and obtain unique fragments
  sort -k9,9 |\
  uniq -u -f8 |\
  # Remove $5 (base pair overlap)
  awk 'OFS="\t" {print $1, $2, $3, $4, $6, $7, $8, $9, $10, $11}' - > tmp_${id}_uniq.bed
  echo "temp_uniq:" tmp_${id}_uniq.bed

  # Intersect fragments with 5Mb bins and obtain duplicated fragments at bin boundaries
  bedtools intersect -a $bins5mb -b $f -wo |\
  # Moved base overlap ($11) to $5 to avoid issues with uniq syntax
  awk 'OFS="\t" {ov=$11; print $1, $2, $3, $4, ov, $5, $6, $7, $8, $9, $10}' - |\
  sort -k9,9 |\
  # Obtain duplicate fragment based QNAME ($9)
  uniq -D -f8 |\
  # Sort by QNAME and base pair overlap ($5)
  sort -k9,9 -k5,5nr |\
  # Compare base pair overlap of duplicate fragment & filter out fragment with smallest overlap
  # If tie, will use the first bin
  awk 'BEGIN {
    FS=OFS="\t"; prev=""; maxval=0} {
        if ($9==prev) {
                if ($5>maxval) {maxval=$5; line=$0}
        } else {
                if (prev!="") print line;
                        maxval=$5; line=$0; prev=$9
        }
  }
  END {if (prev!="") print line}' - |\
  # Remove $5 (base pair overlap)
  awk 'OFS="\t" {print $1, $2, $3, $4, $6, $7, $8, $9, $10, $11}' - > tmp_${id}_dupl.bed
  echo "temp_dupl:" tmp_${id}_dupl.bed

  # Combine unique and de-duplicated fragments
  cat tmp_${id}_uniq.bed tmp_${id}_dupl.bed > $frags_dir/${id}_frags_5mb.bed

  # Delete temporary files
  rm tmp_${id}_uniq.bed tmp_${id}_dupl.bed
  echo done with $f

done
echo " Intersecting fragments with 5Mb bins and resolving duplicates......Finished "



#######################################################################
########### Fragments with size and GC content information ############
# Define output directories
out_dir="${PROJECT_DIR}/05-filter_bedpe_5mb"
mkdir -p "$out_dir" #Create directory if they do not exist


# Iterate through BEDPE files
for f in $(find $bedpe_dir -maxdepth 1 -iname "*.bedpe" -type f)
do
  # Obtain sample name
  id=$(basename -a -s .bedpe $f)
  # Extract fragments from BEDPE file from filtered/binned fragments
  # Appends GC content to 10th column
  awk '
    BEGIN {
      # Set the input and output field separators to tab
      FS=OFS="\t"
    }
    FNR==NR {
      # For the first, create an array with QNAME as the key and GC content as the value
      arr[$8]=$10
      # Skip to the next record without executing the rest of the code
      next
    }
    ($7 in arr) {
      # For the second file (BEDPE), if the QNAME exists in the array,
      # print the current line and append the GC content from the array
      print $0, arr[$7]
    }
  ' $frags_dir/${id}_frags_5mb.bed $f > $out_dir/${id}_filtered.bedpe
done
echo " Extracting filtered/intersected fragments......Finished " 

#######################################################################
################# final fragment with size and GC contents ############

# Define input/output directories
outdir="${PROJECT_DIR}/05-frags_gc"
plotdir="${PROJECT_DIR}/05-gcbias_plots" #GC bias plots
statdir="${PROJECT_DIR}/05-frags_gc_stats" #Filtering statistics

mkdir -p "$outdir" "$plotdir" "$statdir" #Create directories if needed

# Define the path to the R script
R_SCRIPT="${PROJECT_DIR}/05-frags_gc.R"

# Run R script

Rscript "$R_SCRIPT" \
--motifdir $frags_dir \
--outdir $outdir \
--plotdir $plotdir \
--statdir $statdir
echo " performing GC correction and summrazing fragments "

#######################################################################
########################### Combining samples #########################
outdir="${PROJECT_DIR}/05-combine_bin"
mkdir -p "$outdir"  #Create directories if needed

# Define the path to the R script
R_SCRIPT="${PROJECT_DIR}/05-combine_bin.R"

# Initialize R script
Rscript "$R_SCRIPT"  \
--fragdir $fragdir \
--outdir $outdir
echo " Combining data......Finished! "































