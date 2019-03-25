#!/bin/bash
#SBATCH -c 4
#SBATCH -t 1-20:0:0
#SBATCH --mem-per-cpu=16G
#SBATCH -p medium
#SBATCH -o /home/am717/slurm_runs/subbams2.%a.out
#SBATCH -e /home/am717/slurm_runs/subbams2.%a.err
#
# ______________________________________________________________________________________________________________
#
# sbatch --array=1-18 mln_subbamming.sh file.txt /dir/for/outputs/
#
# PARAMETERS:
# 
# 1. file.txt with full pathes for sam files to chop, line by line
# 2. /dir/for/outputs/ is any directory, where all outputs will be placed 
#
# EXAMPLE [for 18 NEB, SMARTseq10ng and SMARTseq100pg data]:
# sbatch --array=1-18 /home/am717/scripts/mln_subbamming.sh /home/am717/sams.txt /n/scratch2/am717/kidney/data/base10_mln_bams
#

module load gcc python R/3.4.1 samtools/1.3.1 star

DIR=$2
mkdir -p $DIR

i=${SLURM_ARRAY_TASK_ID}

##
## 0. From i-th sam from list of files create all possible mln files
##

inputsam=`head -$i $1 | tail -1`
samsize=`samtools view -c $inputsam`

MLNS=`for (( j=$(($samsize / 10000000)); j>=1; j-- )); do echo $(( $j * 10 )); done`

echo "[0] Will create ${MLNS[*]} -mlns files from $inputsam ($samsize reads)."


for mln in $MLNS; do

  declare -a STATS=()
  Nstats=0

  for trial in {1..5}; do

    ###
    ### 1. SAM -> [samtools] -> sampled millions SAM
    ###

    rseed=$(( $RANDOM % 1000 ))
    p=0`bc <<< "scale = 11; ($mln*1000000/$samsize)"`
    s=$( echo $rseed + $p | bc )

    samname=$( basename $inputsam )
    base=${samname%.*}".mln"$mln"_trial"$trial
    ALDIR=$DIR/alignments
    mkdir $ALDIR
    
    CMD11="samtools view -s $s -b -o $ALDIR/${base}.bam $inputsam"
    echo "CMD:    " $CMD11
    $CMD11 

    bam=$ALDIR/$base'.sorted.bam'
    CMD12="samtools sort -O bam -o $bam -@ 4 $ALDIR/${base}.bam"
    echo "CMD:    " $CMD12
    $CMD12

    echo "[1] MLN $mln No $trial : $bam created."

    ##
    ## 2. N mln bam -> [allele counter] -> stat.txt
    ##

    ref129S1="/n/scratch2/am717/references/129S1_pseudo/129S1_pseudo.fa"
    VCF="/n/scratch2/sv111/ASE/mgp.v5.merged.snps_all.dbSNP142_subset_hetero_F1.vcf.gz"

    STDIR=$DIR/stat_allelecouner/input_data
    mkdir -p $STDIR
    stat=$STDIR/$base".stat_0.txt"

    CMD2="python /home/am717/scripts/allelecounter.py --vcf $VCF --sample F1 --bam $bam --ref $ref129S1 --min_cov 0 --min_baseq 2 --min_mapq 10 --o $stat"
    echo "CMD:    " $CMD2 
    $CMD2
    
    STATS+=${base}
    Nstats=$(( Nstats + 1 ))

    echo "[2.$Nstats] Allelecount $stat created."

  done

  ##
  ## 3. stat.txt -> [counts to SNPs + SNPs to genes] -> gene counts tabs
  ##

  echo ${STATS%?}

  CMD_toSNP="Rscript --vanilla /home/am717/scripts/counts_to_SNPs_extended2.R -p $DIR/stat_allelecouner/ -s GRCm38 -e /n/scratch2/sv111/ASE/ -m /n/scratch2/am717/references/F1_pseudo/snp_F1_info_exons.txt"
  CMD_toGenes="Rscript /home/am717/scripts/SNPs_to_genes_extended2.2.R -p $DIR/stat_allelecouner/ -k $Nstats"

  CMD31="$CMD_toSNP -r MLN$mln"_"$base -n $STATS"
  CMD32="$CMD_toGenes -n MLN$mln"_"$base"
  echo "CMD:    " $CMD31
  echo "CMD:    " $CMD32
  $CMD31
  $CMD32

  echo "[3] DataFrame MLN$mln"_"$base at $DIR/stat_allelecouner/ created from $Nstats : ${STATS[*]}."


done




