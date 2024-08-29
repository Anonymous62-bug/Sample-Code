#!/bin/bash

#Get BAM file statistics
exec >> scripts/rseqc.log
exec 2>&1

module purge
ml SAMtools/1.11-GCC-10.2.0
ml RSeQC/3.0.0-foss-2018b-Python-3.6.6
ml R/4.0.0-foss-2019b-fh1

chromsize=/path/to/STAR2Index/chrNameLength.txt
refseq=/path/to/RefSeq.bed
rseq=/path/to/RSeQC/bin


fastq_files=$(ls fastq/*_R1*.fastq.gz)
for file in ${fastq_files[@]};do
    sample=$(echo $(basename $file) | cut -d"_" -f1)
    cd $sample
    inputbam=STAR/Aligned.sortedByCoord.out.bam
    samtools index $inputbam

    /usr/bin/sbatch -N1 -n8 --mem=16G -J BM_${sample} --output=bamstats.report.txt --mail-type=FAIL --mail-user=hsriniva@fredhutch.org --wrap="python $rseq/bam_stat.py  -i $inputbam"
    /usr/bin/sbatch -N1 -n8 --mem=16G -J ID_${sample} --output=log/innerdist.log --mail-type=FAIL --mail-user=hsriniva@fredhutch.org --wrap="python $rseq/inner_distance.py -i $inputbam -o output -r $refseq"
    /usr/bin/sbatch -N1 -n8 --mem=16G -J RD_${sample} --output=readDist.report.txt --mail-type=FAIL --mail-user=hsriniva@fredhutch.org --wrap="python $rseq/read_distribution.py  -i $inputbam -r $refseq"
    cd ..
done
