
#########################################################
# 							#
# 		ASE Replicates Project			#
#         https://github.com/gimelbrantlab/ASE		#
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
# 
# DEPEND:
# bash(v4.2.46)
# gcc(v6.2.0)
# python(v3.6.0)
# vcftools(v0.1.13)
# STAR(v2.5.4a)
#
# Please, use --help|-h for help. 
# For more information, please, visit:
#       https://github.com/gimelbrantlab/ASE
#

import argparse
import os
import tempfile
import subprocess
import gzip

def GATK_SnpExon_Caller(ref, vcf, ofile):
    '''returns the GATK command for VCF -> SnpExon VCF'''
    formulae = ["java -Xmx4g -jar GenomeAnalysisTK.jar -T SelectVariants -R", "-V", "-o", "-selectType SNP -restrictAllelesTo BIALLELIC"]
    cmd = " ".join([formulae[0], ref, formulae[1], vcf, formulae[2], ofile, formulae[3]]
    subprocess.check_output(cmd, shell=True)
    print(cmd) 
    return
def GATK_SepSnpExon_Caller(ref, vcf, ofile, name):
    '''returns the GATK command for joint VCF -> separate SnpExon VCF'''
    formulae = ["java -Xmx4g -jar GenomeAnalysisTK.jar -T SelectVariants -R", "-V", "-o", "-sn", "-selectType SNP -restrictAllelesTo BIALLELIC"]
    cmd = " ".join([formulae[0], ref, formulae[1], vcf, formulae[2], ofile, formulae[3], name, formulae[4]]
    subprocess.check_output(cmd, shell=True)
    print(cmd)
    return

def gzip_tabix_VCF(vcf):
    '''bgzip+tabix'''
    cmd_bgzip = " ".join(["bgzip -c", vcf, ">", vcf + '.gz'])
    subprocess.check_output(cmd_bgzip, shell=True)
    print(cmd_bgzip)
    cmd_tabix = " ".join(["tabix -p vcf", vcf + '.gz']),
    subprocess.check_output(cmd_tabix, shell=True)
    print(cmd_tabix)
    return

def main():
    # ARGUMENTS HANDLING: 
    parser = argparse.ArgumentParser()
    parser.add_argument("--PSEUDOREF", required=False, default="True", metavar="[True/False]", help="Pseudogenomes creation is needed? Default: True")
    parser.add_argument("--F1VCF", required=False, default="True", metavar="[True/False]", help="F1 vcf creation is needed? Default: True")
    parser.add_argument("--pseudo_dir", required=True, metavar="[/path/to/dir]", help="Path to directory with subdirectories --name_mat and --name_pat, that contain pseudogenome fasta files.")
    parser.add_argument("--f1_dir", required=True, metavar="[/path/to/dir]", help="Path to directory with F1 files.")
    parser.add_argument("--ref", required=True, metavar="[/path/to/file]", help="Path to reference fasta file.")
    parser.add_argument("--gtf", required=True, metavar="[/path/to/file]", help="Path to reference gtf file.")
    parser.add_argument("--name_mat", required=False, default="mat", help="Name of maternal line/sample (required with --vcf_joint: names should coincide with names in vcf; if not provided: 'mat')")
    parser.add_argument("--name_pat", required=False, default="pat", help="Name of paternal line/sample (required with --vcf_joint: names should coincide with names in vcf; if not provided: 'pat')")
    parser.add_argument("--vcf_mat", required=False, metavar="[/path/to/file]", help="Path to maternal vcf file (eigther --vcf_joint or pair --vcf_mat & --vcf_pat required; separate vcf will dominate if both provided)")
    parser.add_argument("--vcf_pat", required=False, metavar="[/path/to/file]", help="Path to paternal vcf file (eigther --vcf_joint or pair --vcf_mat & --vcf_pat required; separate vcf will dominate if both provided)")
    parser.add_argument("--vcf_joint", required=False, metavar="[/path/to/file]", help="Path to joint vcf file (eigther --vcf_joint or pair --vcf_mat & --vcf_pat required; separate vcf will dominate if both provided)")
    args = parser.parse_args()

    # TEST IF VCF IS PROVIDED:
    if (args.vcf_joint is None and (args.vcf_mat is None and args.vcf_pat is None)): 
        msg = "Eigther --vcf_joint or pair --vcf_mat & --vcf_pat required."
        raise argparse.ArgumentTypeError(msg) 
    # TEST IF NAMES ARE PROVIDED IN CASE OF JOINT VCF:
    if ((args.vcf_mat is None and args.vcf_pat is None) and (args.name_mat is None or args.name_pat is None)):
        msg = "Required with --vcf_joint: names should coincide with names in vcf."
        raise argparse.ArgumentTypeError(msg)

     # 1. PSEUDOREFERENCE:
     sep_vcf_mat = tempfile.NamedTemporaryFile(delete=False)
     sep_vcf_pat = tempfile.NamedTemporaryFile(delete=False)

     if (args.RSEUDOREF=="True") :
          sep_vcf_mat = tempfile.NamedTemporaryFile(delete=False)
          sep_vcf_pat = tempfile.NamedTemporaryFile(delete=False)

          pseudo_dir_mat = os.path.join(args.pseudo_dir, name_mat)
          pseudo_dir_pat = os.path.join(args.pseudo_dir, name_pat)
          os.makedirs(pseudo_dir_mat, exist_ok=True)
          os.makedirs(pseudo_dir_pat, exist_ok=True)

          # 1.1. Separate SNP VCFs:
          if (args.vcf_mat is None):
              GATK_SepSnpExon_Caller(args.ref, args.vcf_mat, sep_vcf_mat, name_mat)
          else: 
              GATK_SnpExon_Caller(args.ref, args.vcf_mat, sep_vcf_mat)
          if (args.vcf_pat is None):
              GATK_SepSnpExon_Caller(args.ref, args.vcf_pat, sep_vcf_pat, name_pat)
          else:
              GATK_SnpExon_Caller(args.ref, args.vcf_pat, sep_vcf_pat)
          
          gzip_tabix_VCF(sep_vcf_mat)
          gzip_tabix_VCF(sep_vcf_pat)         
          
          os.remove(sep_vcf_mat.name)
          os.remove(sep_vcf_pat.name)


     # 1.2. Pseudoreference:

     # 2. F1 VCF:
     # 2.1 Paired VCF:
     # 2.2 F1 VCF:

if __name__ == "__main__":
    main()



#1. separate.vcf > separate.snp.vcf :
#java -Xmx4g -jar GenomeAnalysisTK.jar -T SelectVariants -R $ref_fa -V $vcf_mat -o "${vcf_mat%.*}".vcf -selectType SNP -restrictAllelesTo BIALLELIC






