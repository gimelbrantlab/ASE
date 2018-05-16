#########################################################
# 							#
# 		ASE Replicates Progect			#
# 		   ..bla-bla-bla..			#
#							#
#   Authors: Mendelevich Asya, Svetlana Vinogradova	#
#							#
#########################################################
#
# DESCRIPTION:
#   Function for full reference preprocessing:
#     - Creating pseudogenomes from reference genomes and vcf file(s)
#     - F1-cross vcf files
#     - Gene-Transcript-Exon Annotations
#Pleasem
# Please, use --help|-h for help.
#

# =====================================================
# 0. ARGUMENTS HANDLING                               ||
# =====================================================
help_message() {
  echo 
  echo "Usage: bash prepare_reference.sh --ref /path/to/reference_genome.fa --gtf /path/to/reference_genome.gtf \ "
  echo "         --name_mat line1 --name_pat line2 --vcf /path/to/joined.vcf.gz \ "
  echo "         --pseudo_dir /path/to/pseudogenomes/output --f1_dir /path/to/f1/files/output --o_dir /path/to/output/dir"
  echo
  echo "  --PSEUDO       if TRUE creates pseudogenomes (default: TRUE; note: if FALSE be shure that subdirectories of --pseudo_dir with pseudogenomes are named the same as maternal and paternal names)"
  echo "  --F1_VCF       if TRUE creates F1 vcf files (default: TRUE; note: if FALSE be shure that F1 vcf file exists with proper name)" 
  echo "  --ref          path to reference_genome.fa file (necessary)"
  echo "  --gtf          path to reference_genome.gtf file (necessary)"
  echo "  --name_mat     name of maternal line/sample (necessary for case of merged vcf: names should coincide with names in vcf; if not provided: 'mat')"
  echo "  --name_pat     name of paternal line/sample (necessary for case of merged vcf: names should coincide with names in vcf; if not provided: 'pat')"
  echo "  --vcf_mat      path to maternal vcf file (eigther --vcf or pair --vcf_mat & --vcf_pat is necessary; separate vcf will dominate if both provided)"
  echo "  --vcf_pat      path to paternal vcf file (eigther --vcf or pair --vcf_mat & --vcf_pat is necessary; separate vcf will dominate if both provided)"
  echo "  --vcf          path to joined vcf file with --name_mat and --name_pat columns (eigther --vcf or pair --vcf_mat & --vcf_pat is necessary; separate vcf will dominate if both provided)"
  echo "  --pseudo_dir   directory for pseudogenomes (necessary)"
  echo "  --f1_dir       directory for F1 vcf files (necessary)"
  echo "  --o_dir        path to output directory (necessary)"
  echo
  echo "Requirements:"
  echo
  echo "  bash(v4.2.46)"
  echo "  gcc(v6.2.0)"
  echo "  python(v3.6.0)"
  echo "  vcftools(v0.1.13)"
  echo "  STAR(v2.5.4a)"
  echo 
  echo "Restrictions: "
  echo
  echo "  * uses '/' as directory path delimiter."
  echo
  echo "If you have more questions, please visit: https://github.com/gimelbrantlab/ASE ."
  echo
  exit 1
}

PSEUDO=TRUE
F1_VCF=TRUE

while [[ $# -gt 0 ]]
do
  case "$1" in 
    --help | -h)
      help_message
      exit
      ;;
    --PSEUDO)
      PSEUDO=$2
      ;;
    --F1_VCF)
      F1_VCF=$2
      ;;
    --pseudo_dir)
      pseudo_dir=$2
      ;;
    --f1_dir)
      f1_dir=$2
      ;;
    --ref)
      ref_fa=$2
      ;;
    --gtf)
      ref_gtf=$2
      ;;
    --name_mat)
      name_mat=$2
      ;;
    --name_pat)
      name_pat=$2
      ;;
    --vcf_mat)
      vcf_mat=$2
      ;;
    --vcf_pat)
      vcf_pat=$2
      ;;
    --vcf)
      vcf_joint=$2
      ;;
    --o_dir)
      o_dir=$2
      ;;
  esac
  shift
  shift
done

raise_variable_absence() {
  echo "ERROR: Not all variables provided. Try --help."
  exit  
}

# Check if everything necessary is provided.
# *** TODO: rewrite (according needs of subscripts) ***

if [[ $ref_fa -eq "" || $ref_gtf -eq "" || $o_dir -eq "" || (($vcf_mat -eq "" || $vcf_pat -eq "") && $vcf_joint -eq "") ]]
then 
  raise_variable_absence
fi

if [[ $name_mat -eq "" ]]
then
  if [[ $vcf_mat -ne "" && $vcf_pat -ne "" ]]
  then
    name_mat="mat"
    echo "WARNING: Maternal name was not provided, set to 'mat'."
  else
    echo "ERROR: Names should be provided in case of joint vcf."
    raise_variable_absence
  fi
fi

if [[ $name_pat -eq "" ]]
then
  if [[ $vcf_mat -ne "" && $vcf_pat -ne "" ]]
  then
    name_mat="pat"
    echo "WARNING: Paternal name was not provided, set to 'pat'"
  else
    echo "ERROR: Names should be provided in case of joint vcf."
    raise_variable_absence
  fi
fi

# =====================================================
# 1. CREATE PSEUDOGENOMES			      ||
# =====================================================

# -----------------------------------------------------
# 1.1. SEPARATE VCFs && SNPs ONLY		      |
# -----------------------------------------------------

if [[ $vcf_mat -eq "" ]]
then
  python vcf_separation.py --vcf_joint --vcf_separate $o_dir"" 
# *** TODO: paths to scripts! ***  
fi 
if [[ $vcf_pat -eq "" ]]
then

fi

# *** TODO: vcfsnp ***

# -----------------------------------------------------
# 1.2. PSEUDOGENOMES                                  |
# -----------------------------------------------------

mkdir $pseudo_dir/$name_mat
mkdir $pseudo_dir/$name_pat
cat $ref_fa | vcf-consensus --vcfsnp_mat  > $pseudo_dir/$name_mat"_pseudo.fa"
cat $ref_fa | vcf-consensus --vcfsnp_pat  > $pseudo_dir/$name_pat"_pseudo.fa"

STAR --runMode genomeGenerate --genomeDir $pseudo_dir/$name_mat --genomeFastaFiles $pseudo_dir/$name_mat"_pseudo.fa" --sjdbGTFfile $ref_gtf
STAR --runMode genomeGenerate --genomeDir $pseudo_dir/$name_pat --genomeFastaFiles $pseudo_dir/$name_pat"_pseudo.fa" --sjdbGTFfile $ref_gtf

# =====================================================
# 1. CREATE F1 VCF FILES                              ||
# =====================================================

# -----------------------------------------------------
# 2.1. F1 VCF, F1 EXON VCF, SNP POSITIONS	      |
# -----------------------------------------------------




