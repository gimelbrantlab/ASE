
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
# NOTE: the order of chromosomes in each file should be the same (add test)!
# FOR MICE ONLY
#
# DEPEND:
# bash(v4.2.46)
# gcc(v6.2.0)
# python(v3.6.0)
# java(v1.8.0_112)
# gatk(v4.0.0.0)
# bcftools(v1.3.1)
# htslib(v1.3.2)
# cufflinks(v2.2.1)
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
import re

def GATK_SelectVariants(r, v, o, g=None, n=None, b=False):
    '''the GATK command for VCF processings; gene - vcf - ofile - gtf - name - biallelic'''
    cmd = "gatk SelectVariants -select-type SNP"
    flags = {'-R':r, '-V':v, '-O':o}
    if (g):
        # tmp exons bed file:
        exon_bed = tempfile.NamedTemporaryFile(delete=False, suffix=".bed")
        cmd_exon = " ".join("grep -w 'exon'", g, "| grep '^[0-9XY]' | awk 'BEGIN{FS=OFS="    "}; {print $1,$4-1,$5}' >", exon_bed.name)
        print(cmd_exon)
        subprocess.check_output(cmd_exon, shell=True)
        flags['-L'] = exon_bed.name
    if (n):
        flags['-sn'] = n
    if (b):
        flags['-restrict-alleles-to'] = "BIALLELIC"
    
    flags_str = ''
    for f in flags:
        if isinstance(flags[f], str):
            flags_str += f + ' ' + flags[f] + ' '
        else:
            for item in flags[f]:
                flags_str += f + ' ' + item + ' '
    cmd = ' '.join([cmd, flags_str])
    print(cmd)
    subprocess.check_output(cmd, shell=True)

    # clear up!:
    if (g): 
        os.remove(exon_bed.name)
    return

def JointToPair_VCF(ref, vcf, ofile, name_mat, name_pat):
    GATK_SelectVariants(r=ref, v=vcf, o=ofile, n=[name_mat, name_pat])
    return

def SepToPair_VCF(ref, vcf_mat, vcf_pat, ofile, name_mat, name_pat):
    unfiltered = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
    cmd_merge = " ".join("bcftools merge", vcf_mat, vcf_pat, "-m both -o", unfiltered.name)
    print(cmd_merge)
    subprocess.check_output(cmd_merge, shell=True)
    GATK_SelectVariants(r=ref, v=unfiltered.name, o=ofile)
    os.remove(unfiltered.name)
    return

def PairToF1_VCF(vcf_pair, vcf_f1, name_mat, name_pat):
    print(vcf_pair + "  -->  " + vcf_f1)
    vcf_stream = open(vcf_pair, 'r')
    out_stream = open(vcf_f1, 'w') 
    
    # header:
    ### REWRITE ACCURATE header filtering ###
    row = vcf_stream.readline()
    while (row.startswith("##")):
        out_stream.write(row)
        row = vcf_stream.readline()
    
    colnames = row.replace('#','').strip().split()
    mat_col = colnames.index(name_mat)
    pat_col = colnames.index(name_pat)
    format_col = colnames.index("FORMAT")
    ref_col = colnames.index("REF")
    alt_col = colnames.index("ALT")
    
    colnames_out = '\t'.join(row[ :min(mat_col, pat_col)] + ["F1"] + row[max(mat_col, pat_col)+1: ])
    out_stream.write(colnames_out)
    
    # body:
    for row in vcf_stream:
        row = row.strip().split()
        gt_index = row[format_col].split(":").index("GT")
        gt_mat = row[mat_col].split(":")[gt_index]
        gt_pat = row[pat_col].split(":")[gt_index]
        if (gt_mat[0]==gt_mat[2] and gt_pat[0]==gt_pat[2] and gt_mat[0]!=gt_pat[0]):
            if (gt_mat[0] == '0'):
                ref_allele = row[ref_col]
            else: 
                ref_allele = row[alt_col].split(',')[gt_mat[0]-1]
            if (gt_pat[0] == '0'):
                alt_allele = row[ref_col]
            else:
                alt_allele = row[alt_col].split(',')[gt_pat[0]-1]
            
            row[pat_col].replace(gt_pat, "0|1")
            row[ref_col] = ref_allele
            row[alt_col] = alt_allele
            ### IF THERE ARE OTHER FIELDS IN COLUMN TO BE CHANGED? ###

            colnames_out = '\t'.join(row[ :min(mat_col, pat_col)] + row[pat_col] + row[max(mat_col, pat_col)+1: ])
            out_stream.write(colnames_out)

    vcf_stream.close()
    out_stream.close()
    return

def gzip_tabix_VCF(vcf):
    '''bgzip+tabix'''
    cmd_bgzip = " ".join(["bgzip -c", vcf, ">", vcf + '.gz'])
    print(cmd_bgzip)
    subprocess.check_output(cmd_bgzip, shell=True)
    cmd_tabix = " ".join(["tabix -p vcf", vcf + '.gz']),
    print(cmd_tabix)
    subprocess.check_output(cmd_tabix, shell=True)
    return

def vcftools_consensus(ref_fa, vcf, pseudo_fa):
    '''creates pseudogenom inserting corresponding SNPs'''
    cmd = " ".join(["cat", ref_fa, "| vcf-consensus", vcf, ">", pseudo_fa])
    print(cmd)
    subprocess.check_output(cmd, shell=True)
    return

def trascriptome_creation():

    # Transcriptomes with cufflinks:
    #gffread /n/scratch2/am717/references/GRCm38/Mus_musculus.GRCm38.68.gtf -g /n/scratch2/am717/references/129S1_pseudo/129S1_pseudo.fa -w /n/scratch2/am717/references/129S1_pseudo_tr/129S1_pseudo_trs.fa
    #gffread /n/scratch2/am717/references/GRCm38/Mus_musculus.GRCm38.68.gtf -g /n/scratch2/am717/references/CAST_pseudo/CAST_pseudo.fa -w /n/scratch2/am717/references/CAST_pseudo_tr/CAST_pseudo_trs.fa
    # Merge transcriptome F1:
    #sed 's/ gene.*$/_129S1/' /n/scratch2/am717/references/129S1_pseudo_tr/129S1_pseudo_trs.fa > /n/scratch2/am717/references/F1_pseudo_tr/129S1_pseudo_trs_marked.fa
    #sed 's/ gene.*$/_CAST/' /n/scratch2/am717/references/CAST_pseudo_tr/CAST_pseudo_trs.fa > /n/scratch2/am717/references/F1_pseudo_tr/CAST_pseudo_trs_marked.fa
    #cat /n/scratch2/am717/references/F1_pseudo_tr/CAST_pseudo_trs_marked.fa /n/scratch2/am717/references/F1_pseudo_tr/129S1_pseudo_trs_marked.fa > /n/scratch2/am717/references/F1_pseudo_tr/F1_pseudo_trs_marked.fa

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
    # SET NAMES:
    if (args.name_mat is None):
        name_mat = "mat"
    else: 
        name_mat = args.name_mat
    if (args.name_pat is None):
        name_pat = "pat"
    else:
        name_pat = args.name_pat

     # 1. PSEUDOREFERENCE:
    if (args.PSEUDOREF=="True"):
        sep_vcf_mat = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
        sep_vcf_pat = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")

        pseudo_dir_mat = os.path.join(args.pseudo_dir, name_mat)
        pseudo_dir_pat = os.path.join(args.pseudo_dir, name_pat)
        os.makedirs(pseudo_dir_mat, exist_ok=True)
        os.makedirs(pseudo_dir_pat, exist_ok=True)

        # 1.1. Separate SNP VCFs:
        if (args.vcf_mat is None):
            GATK_SelectVariants(r=args.ref, v=args.vcf_joint, o=sep_vcf_mat.name, n=name_mat, b=True)
        else: 
            GATK_SelectVariants(r=args.ref, v=args.vcf_mat, o=sep_vcf_mat.name, b=True)
        if (args.vcf_pat is None):
            GATK_SelectVariants(r=args.ref, v=args.vcf_joint, o=sep_vcf_pat.name, n=name_pat, b=True)
        else:
            GATK_SelectVariants(r=args.ref, v=args.vcf_pat, o=sep_vcf_pat.name, b=True)
          
        gzip_tabix_VCF(sep_vcf_mat.name)
        gzip_tabix_VCF(sep_vcf_pat.name)         
          
        # 1.2. Pseudoreference:
        vcftools_consensus(args.ref, sep_vcf_mat.name+'.gz', os.path.join(pseudo_dir_mat, name_mat + "_pseudo.fa"))
        vcftools_consensus(args.ref, sep_vcf_pat.name+'.gz', os.path.join(pseudo_dir_pat, name_pat + "_pseudo.fa"))
         
        os.remove(sep_vcf_mat.name); os.remove(sep_vcf_mat.name+'.gz'); os.remove(sep_vcf_mat.name+'.gz.tbi')
        os.remove(sep_vcf_pat.name); os.remove(sep_vcf_pat.name+'.gz'); os.remove(sep_vcf_pat.name+'.gz.tbi')

    # 2. F1 VCF:
    if (args.F1VCF=="True") :
        os.makedirs(f1_dir, exist_ok=True)
         
        pair_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
         
        # 2.1 Paired VCF:
        if (args.vcf_mat is None or args.vcf_pat is None):
            JointToPair_VCF(args.ref, args.vcf_joint, pair_vcf.name, name_mat, name_pat)
        else: 
            SepToPair_VCF(args.ref, args.vcf_mat, args.vcf_pat, pair_vcf.name, name_mat, name_pat)
         
        # 2.2 F1 VCF:
        vcf_f1 = os.path.join(f1_dir, "_".join("F1", name_mat, name_pat)+'.vcf')
        PairToF1_VCF(pair_vcf.name, vcf_f1, name_mat, name_pat)
        gzip_tabix_VCF(vcf_f1)

        # 2.3 Exon F1 VCF:
        vcf_f1exon = os.path.join(f1_dir, "_".join("F1", name_mat, name_pat, "exon")+'.vcf')
        GATK_SelectVariants(r=args.ref, v=vcf_f1, g=args.gtf, o=vcf_f1exon)
        gzip_tabix_VCF(vcf_f1exon)

if __name__ == "__main__":
    main()

