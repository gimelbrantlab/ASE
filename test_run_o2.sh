#!/bin/bash
#SBATCH -c 8
#SBATCH -t 1-10:0:0
#SBATCH --mem-per-cpu=8G
#SBATCH -p medium
#SBATCH -J test_run_o2
#SBATCH -o /home/am717/slurm_runs/%J.%a.out
#SBATCH -e /home/am717/slurm_runs/%J.%a.err

module load gcc python/3.6.0 java gatk/4.0.0.0 bcftools star samtools picard htslib

# Reference preparation:
python3 /home/am717/ASE/prepare_reference.py --f1_dir /n/scratch2/am717/ASE_test/F1_fromjoint/ --pseudo_dir /n/scratch2/am717/ASE_test/pseudoref_fromjoint/ --ref /n/scratch2/sv111/ASE/GRCm38_68.fa --gtf /n/scratch2/am717/ASE_test/data/Mus_musculus.GRCm38.68_charchrorder.gtf --name_mat 129S1_SvImJ --name_pat CAST_EiJ --vcf_joint /n/scratch2/sv111/ASE/mgp.v5.merged.snps_all.dbSNP142.vcf.gz

# RNAseq files preparation:
#python3 /home/am717/ASE/prepare_RNAseq.py --method CASTEL --fastq_dir /n/sratch2/sv111/data/kidney/data/ --o_dir /n/scratch2/am717/ASE_test/STAR_CASTEL/ --genome /n/scratch2/am717/ASE_test/pseudoref_fromjoint/129S1_SvImJ/*.fa,/n/scratch2/am717/ASE_test/pseudoref_fromjoint/CAST_EiJ/*.fa
