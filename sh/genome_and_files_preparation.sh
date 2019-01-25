#!/bin/bash
#
#BSUB -L /bin/bash             
#BSUB -W 10:00               
#BSUB -q priority       
#BSUB -J starInd                    
#BSUB -o myjob.%J.out            
#BSUB -e myjob.%J.err   
#BSUB -n 6   

#  runPipeline.sh
#  
#
#  Created by Svetlana on 02/26/16
#

# Data is downloaded from ftp://ftp-mouse.sanger.ac.uk/

source /opt/Modules/3.2.10/init/bash
module purge
module load dev/openblas/0.2.14
module load dev/python/2.7.6
module load seq/pysam/0.8.4
module load stats/R/3.3.1
module load dev/java/jdk1.8
module load seq/tabix/0.2.6
export PERL5LIB=/home/sv111/tools/vcftools_0.1.13/perl/
module load seq/samtools/1.3
module load seq/STAR/2.5.2a
export OMP_NUM_THREADS=6

# Download individual SNPs
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz.tbi
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz.tbi

#Download ref fasta and annotation file
wget ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa
wget ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz
gunzip Mus_musculus.GRCm38.68.gtf.gz

## Create pseudogenomes
cat GRCm38_68.fa | /home/sv111/tools/vcftools_0.1.13/bin/vcf-consensus CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz  > CAST_pseudo.fa
cat GRCm38_68.fa | /home/sv111/tools/vcftools_0.1.13/bin/vcf-consensus 129S1_SvImJ.mgp.v5.snps.dbSNP142.vcf.gz  > 129S1_pseudo.fa
mkdir 129S1SvImJ_pseudo/
mkdir CastEij_pseudo/
mv 129S1_pseudo.fa 129S1SvImJ_pseudo/
mv CAST_pseudo.fa CastEij_pseudo/

# Index pseudogenome files with STAR
STAR --runMode genomeGenerate --genomeDir /n/scratch2/v111/SE/CastEij_pseudo --genomeFastaFiles /n/scratch2/sv111/ASE/CastEij_pseudo/CAST_pseudo.fa --sjdbGTFfile /n/scratch2/sv111/ASE/Mus_musculus.GRCm38.68.gtf
STAR --runMode genomeGenerate --genomeDir /n/scratch2/sv111/ASE/129S1SvImJ_pseudo --genomeFastaFiles /n/scratch2/sv111/ASE/129S1SvImJ_pseudo/129S1_pseudo.fa --sjdbGTFfile /n/scratch2/sv111/ASE/Mus_musculus.GRCm38.68.gtf

## Download all SNPs 
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz.tbi

# Convert all SNPs file to SNP file with two strains
java -Xmx4g -jar GenomeAnalysisTK.jar -T SelectVariants -R GRCm38_68.fa -V mgp.v5.merged.snps_all.dbSNP142.vcf.gz -o mgp.v5.merged.snps_all.dbSNP142_subset.vcf -sn CAST_EiJ -sn 129S1_SvImJ -selectType SNP -restrictAllelesTo BIALLELIC

## Gzip and index resulting vcf file
/home/sv111/tools/htslib-1.3.2/bgzip -c mgp.v5.merged.snps_all.dbSNP142_subset.vcf > mgp.v5.merged.snps_all.dbSNP142_subset.vcf.gz
/home/sv111/tools/htslib-1.3.2/tabix -p vcf mgp.v5.merged.snps_all.dbSNP142_subset.vcf.gz

# Filter vcf file to include only heterozygous positions and create F1 column with 0/1 tag
python ./modVcf.py --vcf mgp.v5.merged.snps_all.dbSNP142_subset.vcf.gz --o mgp.v5.merged.snps_all.dbSNP142_subset_hetero_F1.vcf
/home/sv111/tools/htslib-1.3.2/bgzip -c mgp.v5.merged.snps_all.dbSNP142_subset_hetero_F1.vcf > mgp.v5.merged.snps_all.dbSNP142_subset_hetero_F1.vcf.gz
/home/sv111/tools/htslib-1.3.2/tabix -p vcf mgp.v5.merged.snps_all.dbSNP142_subset_hetero_F1.vcf.gz
