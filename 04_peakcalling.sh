#!/bin/bash

Input_path="path to Input files"
Treat_path="path to ChIP files"
cd ${Input_path}
Input_files=`ls *bam`
cd ${Treat_path}
Treat_files=`ls *bam`

for treat in $Treat_files;do
#####	Input="Input"
	for input in ${Input_files};do
		if [ "${input%-*}" = "${treat%-*}" ];then
 			sicer -t ${Treat_path}/${treat} -c ${Input_path}/${input} -s hg38
			awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$4,$5,$6}' ${Treat_path}/${treat%.*}.bed  > ${Treat_path}/${treat%.*}.bed.chr
			mv -f ${Treat_path}/${treat%.*}.bed.chr ${Treat_path}/${treat%.*}.bed
			awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$4,$5,$6}' ${Input_path}/${input%.*}.bed  > ${Input_path}/${input%.*}.bed.chr
			mv -f ${Input_path}/${input%.*}.bed.chr ${Input_path}/${input%.*}.bed
			sicer -t ${Treat_path}/${treat%.*}.bed  -c ${Input_path}/${input%.*}.bed -s hg38 -w 200 -f 150 -egf 0.88 -fdr 0.001 -g 400 -cpu 100 -rt 1
	
		fi
	done

done

