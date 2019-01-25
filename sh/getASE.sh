#!/bin/bash
#
#BSUB -L /bin/bash
#BSUB -W 200:00
#BSUB -q priority
#BSUB -J starInd
#BSUB -o myjob.%J.out
#BSUB -e myjob.%J.err
#BSUB -n 4
#BSUB -R "rusage[mem=85000]"

# ------------------------------------------------------------------
# Svetlana Vinogradova 
#          Allelic Imbalance Pipeline
#
#          This script runs STAR mapping,
#          merges files with maternal and paternal mappings,
#          runs mpileup to count number of counts per SNP per allele
#
# ------------------------------------------------------------------

source /opt/Modules/3.2.10/init/bash
module load seq/STAR/2.5.2a
module load seq/samtools/1.3
module load seq/fastqc/0.11.3
module purge
module load dev/openblas/0.2.14
module load dev/python/2.7.6
module load seq/pysam/0.8.4
module load stats/R/3.3.1
module load seq/tabix/0.2.6
module load seq/samtools/1.3

GENOMEDIRPATH=/n/scratch2/sv111/ASE/
GTF=/n/scratch2/sv111/ASE/Mus_musculus.GRCm38.68.gtf
PATH=/n/scratch2/sv111/data/Abelson/

clones=( "4.11" "H8" "Abl1" "Abl2" "Abl3" )
for name in ${clones[@]}
do
# run STAR: mapping to paternal genome (Cast)
	STAR --genomeDir $GENOMEDIRPATH"CastEij_pseudo" \
    --runThreadN 4 \
    --outFilterMultimapNmax 1 \
    --outSAMtype SAM SortedByCoordinate \
    --readFilesIn $PATH"data/"$name".fastq" \
    --outFileNamePrefix $PATH"results/"$name"ToCast" \
    --outSAMattrRGline ID:pat \
    --sjdbGTFfile GTF
# run STAR: mapping to maternal genome (129)        
    STAR --genomeDir $GENOMEDIRPATH"129S1SvImJ_pseudo" \
    --runThreadN 4 \
    --outFilterMultimapNmax 1 \
    --outSAMtype SAM SortedByCoordinate \
	--readFilesIn $PATH"data/"$name".fastq" \
    --outFileNamePrefix $PATH"results/"$name"To129S1" \ 
    --outSAMattrRGline ID:mat \
    --sjdbGTFfile GTF
# run python to merge two sam files
	python $GENOMEDIRPATH"alleleseq_merge.py" --pat_sam $PATH"results/"$name"ToCastAligned.sortedByCoord.out.sam" --mat_sam $PATH"results/"$name"To129S1Aligned.sortedByCoord.out.sam" --o $PATH"results/"$name"_merged.sam"
# run samtools to convert for further use
	samtools view -Sb $PATH"results/"$name"_merged.sam" > $PATH"results/"$name"_merged.bam"
	samtools sort $PATH"results/"$name"_merged.bam" -T $name"_mergedUpd" > $PATH"results/"$name"_merged.sorted.bam"  
# run ASE counter
	python $GENOMEDIRPATH"allelecounter.py" --vcf $GENOMEDIRPATH"mgp.v5.merged.snps_all.dbSNP142_subset_hetero_F1.vcf.gz" --sample F1 --bam $PATH"results/"$name"_merged.sorted.bam" --ref $GENOMEDIRPATH"129S1SvImJ_pseudo/129S1_pseudo.fa" --min_cov 0 --min_baseq 2 --min_mapq 10 --o $PATH"results/"$name"_stat.txt"    
done
