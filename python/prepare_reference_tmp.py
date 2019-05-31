
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
    '''
    The GATK command for VCF processings, selects SNPs from all variants
    (optionally, restricts on exons only, or takes only biallelic variants)
    Input:  
    r -- genome.fa (should be indexed (samtools faidx) and have dictionary file (Picard CreateSequenceDictionary))
    v -- vcf file (with one column and with only biallelic variants for particular organism (not a set) remains)
    o -- ofile to place output 
    g -- (optional) gtf file for exon positions annotation
    n -- (optional) name of the column to chop from mixed vcf 
    b -- (optional) restriction to biallelic variants flag
    Output: 
    returns nothing
    creates vcf file (name defined via o option) with selected variants
    '''

    cmd = "gatk SelectVariants -select-type SNP "
    flags = {'-R':r, '-V':v, '-O':o}
    if (g):
        # tmp exons bed file:
        exon_bed = tempfile.NamedTemporaryFile(delete=False, suffix=".bed")
        cmd_exon = " ".join(["grep -w 'exon'", g, "| grep '^[0-9XY]' | awk 'BEGIN{FS=OFS=", '"\t"', "}; {print $1,$4-1,$5}' >", exon_bed.name])
        print(cmd_exon)
        subprocess.check_output(cmd_exon, shell=True)
        flags['-L'] = exon_bed.name
    if (n):
        flags['-sn'] = n
    if (b):
        flags['--restrict-alleles-to'] = "BIALLELIC"

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

def ParentalSeparation_VCF(v, o_1, o_2, ind_name):
    '''
    Takes individual fazed vcf and 
    Input:
    v  -- path to individual vcf
    o1 -- path to output file for first allele vcf
    o2 -- path to output file for second allele vcf
    ind_name -- name of individual, should coinside with column name in vcf
    Output:
    Two vcf files for two alleles of the organism, reference is reference, haplotype 1|1, for each allele
    '''
    vcf_stream = open(v, 'r')
    out1_stream = open(o_1, 'w')
    out2_stream = open(o_2, 'w')

    # Read header:
    row = vcf_stream.readline()
    while (row.startswith("##")):
        out1_stream.write(row)
        out2_stream.write(row)
        row = vcf_stream.readline()

    # Column Names: 
    colnames = row.replace('#','').strip().split()
    ind_col = colnames.index(ind_name)
    format_col = colnames.index("FORMAT")
    ref_col = colnames.index("REF")
    alt_col = colnames.index("ALT")

    colnames[ind_col] = ind_name + ".mat"
    out1_stream.write('#' + '\t'.join(colnames) + '\n')
    colnames[ind_col] = ind_name + ".pat"
    out2_stream.write('#' + '\t'.join(colnames) + '\n')

    # Row by row:
    for row in vcf_stream:
        row = row.strip().split()
        gt_index = row[format_col].split(":").index("GT")
        gt_ind_list = row[ind_col].split(":")
        gt_ind = row[ind_col].split(":")[gt_index]

        if (len(gt_ind)==1 and gt_ind!='.'):
            
            if (gt_ind == '0'):
                gt_ind_12 = '0|0'
                alt_allele_12 = row[ref_col]
            else:
                gt_ind_12 = '1|1'
                alt_allele_12 = row[alt_col].split(',')[int(gt_ind)-1]

            gt_ind_list[gt_index] = gt_ind_12
            row[ind_col] = ":".join(gt_ind_list)
            row[alt_col] = alt_allele_12
            cols_out12 = '\t'.join(row)
            out1_stream.write(cols_out12 + '\n')
            out2_stream.write(cols_out12 + '\n')

        elif (gt_ind[0]!='.' and gt_ind[2]!='.'):

            if (gt_ind[0] == '0'):
                gt_ind_1 = '0|0'
                alt_allele_1 = row[ref_col]
            else:
                gt_ind_1 = '1|1'
                alt_allele_1 = row[alt_col].split(',')[int(gt_ind[0])-1]

            if (gt_ind[2] == '0'):
                gt_ind_2 = '0|0'
                alt_allele_2 = row[ref_col]
            else:
                gt_ind_2 = '1|1'
                alt_allele_2 = row[alt_col].split(',')[int(gt_ind[2])-1]

            gt_ind_list[gt_index] = gt_ind_1
            row[ind_col] = ":".join(gt_ind_list)
            row[alt_col] = alt_allele_1
            cols_out1 = '\t'.join(row)
            out1_stream.write(cols_out1 + '\n')

            gt_ind_list[gt_index] = gt_ind_2
            row[ind_col] = ":".join(gt_ind_list)
            row[alt_col] = alt_allele_2
            cols_out2 = '\t'.join(row)
            out2_stream.write(cols_out2 + '\n')

        ### IF THERE ARE OTHER FIELDS IN COLUMN TO BE CHANGED? ###

    vcf_stream.close()
    out1_stream.close()
    out2_stream.close()
    return

def takeBiallelic(vcf):


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
    '''creates pseudogenom inserting corresponding SNPs (vcftools vcf-consensus)'''
    # it is VERY SLOW, maybe replace with something else?
    #cmd = " ".join(["<", ref_fa, "vcf-consensus", vcf, ">", pseudo_fa])
    cmd = " ".join(["cat", ref_fa, "| vcf-consensus", vcf, ">", pseudo_fa])
    print(cmd)
    subprocess.check_output(cmd, shell=True)
    return


def main():
    # ARGUMENTS HANDLING: 
    parser = argparse.ArgumentParser()
    parser.add_argument("--PSEUDOREF", required=False, default="False", metavar="[True/False]", help="Pseudogenomes creation is needed? Default: False")
    parser.add_argument("--HETVCF", required=False, default="False", metavar="[True/False]", help="Het vcf creation is needed? Default: False")
    parser.add_argument("--pseudoref_dir", required=False, metavar="[/path/to/dir]", help="Path to directory with subdirectories --name_mat and --name_pat, that contain pseudogenome fasta files.")
    parser.add_argument("--vcf_dir", required=False, metavar="[/path/to/dir]", help="Path to directory with vcf files; required if --HETVCF True, or allelevcfs not provided")
    parser.add_argument("--ref", required=True, metavar="[/path/to/file]", help="Path to reference fasta file.")
    parser.add_argument("--gtf", required=False, metavar="[/path/to/file]", help="Path to reference gtf file.")
    parser.add_argument("--name_ind", required=False, help="Name of individuum (required with --vcf_joind and --vcf_ind: name should coincide with name in vcf; please avoid giving 'ref' ar 'alt' names, they have special meanings)")
    parser.add_argument("--name_mat", required=False, help="Name of maternal line/sample (required with --vcf_joint if aat is needed, and --vcf_mat: names should coincide with names in vcf; please avoid giving 'ref' ar 'alt' names, they have special meanings)")
    parser.add_argument("--name_pat", required=False, help="Name of paternal line/sample (required with --vcf_joint if pat is needed, and --vcf_pat: names should coincide with names in vcf; please avoid giving 'ref' ar 'alt' names, they have special meanings)")
    parser.add_argument("--vcf_mat", required=False, metavar="[/path/to/file]", help="Path to maternal vcf file (if F1 cross, eigther --vcf_joint or pair --vcf_mat & --vcf_pat required; separate vcf will dominate if both provided)")
    parser.add_argument("--vcf_pat", required=False, metavar="[/path/to/file]", help="Path to paternal vcf file (if F1 cross, eigther --vcf_joint or pair --vcf_mat & --vcf_pat required; separate vcf will dominate if both provided)")
    parser.add_argument("--vcf_joint", required=False, metavar="[/path/to/file]", help="Path to joint vcf file for lines (if F1 cross, eigther --vcf_joint or pair --vcf_mat & --vcf_pat required; separate vcf will dominate if both provided)")
    parser.add_argument("--vcf_ind", required=False, metavar="[/path/to/file]", help="Path to individuum vcf file (if not F1 cross, eigther --vcf_joind or --vcf_ind required; separate vcf will dominate if both provided)")
    parser.add_argument("--vcf_joind", required=False, metavar="[/path/to/file]", help="Path to joint vcf file for individuums (if not F1 cross, eigther --vcf_joind or --vcf_ind required; separate vcf will dominate if both provided)") 

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # TEST IF EVERYTHING NECESSARY IS PRESENT and SET NAMES and A MODE |
    # ------------------------------------------------------------------
    # GENERAL PARAMs:

    if (args.ref is None):
         msg = "Required parameter --ref is missing."
         raise argparse.ArgumentTypeError(msg)
    if (args.PSEUDOREF is True and args.pseudoref_dir is None):
         msg = "Required parameter --pseudoref_dir is missing."
         raise argparse.ArgumentTypeError(msg)
    # if (args.HETVCF is True and args.vcf_dir is None):
    #      msg = "Required parameter --vcf_dir is missing."
    #      raise argparse.ArgumentTypeError(msg)
    if (args.vcf_dir is None):
         msg = "Required parameter --vcf_dir is missing."

    # CASES:
    # SEPARATE ALLELE VCFs:
    if (args.vcf_mat is not None or args.vcf_pat is not None):
        if ((args.vcf_mat is not None and args.vcf_mat is not None) and (args.vcf_pat is not None and args.vcf_pat is not None)):
            input_case = "two_alleles"
            name_mat   = args.name_mat
            name_pat   = args.name_pat
        elif (args.vcf_mat is not None and args.vcf_mat is not None):
            input_case = "one_allele"
            vcf_alt    = args.vcf_mat
            name_alt   = args.name_mat
        elif (args.vcf_pat is not None and args.vcf_pat is not None):
            input_case = "one_allele"
            vcf_alt    = args.vcf_pat
            name_alt   = args.name_pat
        else: 
            msg = "Data required: --vcf_xxx should go in a pair with --name_xxx."
            raise argparse.ArgumentTypeError(msg)
    # JOINT VCF FILE WITH ALLELES as columns:
    elif (args.vcf_joint is not None):
        if (args.vcf_mat is not None and args.vcf_pat is not None):
            input_case = "joint_two_alleles"
            name_mat   = args.name_mat
            name_pat   = args.name_pat
        elif (args.vcf_mat is not None):
            input_case = "joint_one_allele"
            name_alt   = args.name_mat
        elif (args.vcf_pat is not None):
            input_case = "joint_one_allele"
            name_alt   = args.name_pat
        else:
            msg = "Data required: at least either --name_mat or --name_pat is needed."
            raise argparse.ArgumentTypeError(msg)
    # INDIVIDUAL or JOINT INDUVUDUAL VCF (not a cross):
    elif (args.vcf_ind is not None or args.vcf_joind is not None):
        if (args.name_ind is not None):
            name_ind   = args.name_ind
            name_mat   = str(args.name_ind) + "_mat"
            name_pat   = str(args.name_ind) + "_pat"
            if (args.vcf_ind is not None):
                input_case = "individ"
            elif (args.vcf_joind is not None):
                input_case = "joint_individs"
        else:
            msg = "Data required: parameter --name_ind is missing."
            raise argparse.ArgumentTypeError(msg)
    # SMTH INCORRECT:
    else:
        msg = "Check the correctness of required parameters for your case. See help."
        parser.print_help()
        raise argparse.ArgumentTypeError(msg)

    # ------------------------------------------------------------------
    # PSEUDOREFERENCE CREATION: if PSEUDOREF set to be True            |
    # ------------------------------------------------------------------
 
    if (args.PSEUDOREF=="True"):

        cmd_mkdir_vcf = "mkdir -p " + args.vcf_dir
        subprocess.check_output(cmd_mkdir_vcf, shell=True)

        # TMP VCFs an DIRECTORIES:
        if (input_case == "two_alleles" or input_case == "joint_two_alleles" or input_case == "individ" or input_case == "joint_individs"):
            # Output vcfs:
            sep_vcf_mat = os.path.join(args.vcf_dir, name_mat + ".SNP.biallelic.vcf")
            sep_vcf_pat = os.path.join(args.vcf_dir, name_pat + ".SNP.biallelic.vcf")
            # Creation of directories for pseudoref:
            pseudo_dir_mat = os.path.join(args.pseudoref_dir, name_mat)
            pseudo_dir_pat = os.path.join(args.pseudoref_dir, name_pat)
            cmd_mkdir_mat = "mkdir -p " + pseudo_dir_mat
            subprocess.check_output(cmd_mkdir_mat, shell=True)
            cmd_mkdir_pat = "mkdir -p " + pseudo_dir_pat
            subprocess.check_output(cmd_mkdir_pat, shell=True)
            # All cases:
            if (input_case == "two_alleles"):
                GATK_SelectVariants(r=args.ref, v=args.vcf_mat, o=sep_vcf_mat, b=True)
                GATK_SelectVariants(r=args.ref, v=args.vcf_pat, o=sep_vcf_pat, b=True)
            elif (input_case == "joint_two_alleles"):
                GATK_SelectVariants(r=args.ref, v=args.vcf_joint, o=sep_vcf_mat, n=name_mat, b=True)
                GATK_SelectVariants(r=args.ref, v=args.vcf_joint, o=sep_vcf_pat, n=name_pat, b=True)
            elif (input_case == "individ" or input_case == "joint_individs"):
                #sep_vcf_ind = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf")
                sep_vcf_ind = os.path.join(args.vcf_dir, name_ind + ".individual.SNP.vcf")
                if (input_case == "individ"):
                    GATK_SelectVariants(r=args.ref, v=args.vcf_ind, o=sep_vcf_ind, b=False)
                elif (input_case == "joint_individs"):
                    GATK_SelectVariants(r=args.ref, v=args.vcf_joind, o=sep_vcf_ind, n=name_ind, b=False)
                ParentalSeparation_VCF(v=sep_vcf_ind, o_1=sep_vcf_mat, o_2=sep_vcf_pat, ind_name=name_ind)
                #os.remove(sep_vcf_ind.name)
            # Indexing vcfs:
            gzip_tabix_VCF(sep_vcf_mat)
            gzip_tabix_VCF(sep_vcf_pat)
            # Pseudoreference:
            vcftools_consensus(args.ref, sep_vcf_mat + '.gz', os.path.join(pseudo_dir_mat, name_mat + "_pseudo.fa"))
            vcftools_consensus(args.ref, sep_vcf_pat + '.gz', os.path.join(pseudo_dir_pat, name_pat + "_pseudo.fa"))

        elif (input_case == "one_allele" or input_case == "joint_one_allele"):
            # Alt vcfs:
            sep_vcf_alt = os.path.join(args.vcf_dir, name_alt + ".SNP.biallelic.vcf")
            # Creation of directories for pseudoref:
            pseudo_dir_alt = os.path.join(args.pseudoref_dir, name_alt)
            cmd_mkdir_alt = "mkdir -p " + pseudo_dir_pat
            subprocess.check_output(cmd_mkdir_alt, shell=True)
            # All cases:
            if(input_case == "one_allele"):
                GATK_SelectVariants(r=args.ref, v=vcf_alt, o=sep_vcf_alt, b=True)
            elif(input_case == "joint_one_allele"):
                GATK_SelectVariants(r=args.ref, v=args.vcf_joint, n=name_alt, o=sep_vcf_alt, b=True)
            # Indexing vcfs:
            gzip_tabix_VCF(sep_vcf_alt)
            # Pseudoreference:
            vcftools_consensus(args.ref, sep_vcf_alt + '.gz', os.path.join(pseudo_dir_alt, name_alt + "_pseudo.fa"))

        else:
            msg = "Something went wrong in the input data parsing."
            raise argparse.ArgumentTypeError(msg)





if __name__ == "__main__":
    main()

