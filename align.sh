#!/bin/bash

#Align FASTQ reads to reference genome
exec >> "scripts/align.log"
exec 2>&1

module purge
module load STAR/2.7.7a-GCC-10.2.0

fastq_files=$(ls fastq/*_R1*.fastq.gz)

STARref=/path/to/STAR2Index
gtfFile=/path/to/gencode.annotation.gtf

for file in ${fastq_files[@]};do
    sample=$(echo $(basename $file) | cut -d"." -f1 | sed 's#_R1##g')
    r1=${file}
    r2=$( echo ${r1} | sed 's#R1#R2#')
    mkdir -p ${sample}/STAR
    cd ${sample}/STAR

    #START JOB
    job_id=$(/usr/bin/sbatch -N1 -n8 --mem=32G -J align_${sample} --output=../log/%x.log --wrap="time STAR \
    --genomeDir $STARref \
    --readFilesIn ../../${r1} ../../${r2} \
    --readFilesCommand zcat \
    --runThreadN 4 \
    --sjdbGTFfile $gtfFile \
    --outFilterMultimapNmax 20 \
    --alignIntronMax 500000 \
    --alignMatesGapMax 1000000 \
    --sjdbScore 2 \
    --alignSJDBoverhangMin 1 \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --sjdbOverhang 100 \
    --twopassMode Basic \
    --outSAMstrandField intronMotif \
    --outSAMattributes NH HI NM MD AS XS \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMheaderHD @HD VN:1.4 \
    --quantMode TranscriptomeSAM GeneCounts ")
    cd ../../
done
