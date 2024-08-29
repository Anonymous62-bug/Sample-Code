#!/bin/bash

## FASTP program to look into fastq file data quality
exec >> "scripts/fastp.log"
exec 2>&1

module purge
ml fastp/0.20.0-GCC-8.3.0

fastq_files=$(ls fastq/*_R1*.fastq.gz)

for file in ${fastq_files[@]};do
    sample=$(echo $(basename $file) | cut -d"." -f1 | sed 's#_R1##g')
    [[ -d ${sample} ]] || mkdir ${sample}
    r1=${file}
    r2=$( echo ${r1} | sed 's#R1#R2#')

    #Start job
    job_id=$(/usr/bin/sbatch -N1 -n8 --mem=32G -J fastp_${sample} --output=${sample}/%x.report.txt --wrap="time \
    fastp \
    -o fastp \
    --in1 ${r1} \
    --in2 ${r2} \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --overrepresentation_analysis \
    --overrepresentation_sampling 20 \
    --compression 4 \
    --out1 ${sample}/temp_${sample}_R1.fastq.gz \
    --out2 ${sample}/temp_${sample}_R2.fastq.gz \
    --report_title ${sample} \
    --json ${sample}/fastp.json \
    --html ${sample}/fastp.html")
done
