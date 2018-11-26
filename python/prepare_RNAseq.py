
#########################################################
#                                                       #
#               ASE Replicates Project                  #
#         https://github.com/gimelbrantlab/ASE          #
#                                                       #
#   Authors: Mendelevich Asya, Svetlana Vinogradova     #
#                                                       #
#########################################################
#
# DESCRIPTION:
#   Function for RNAseq fastq files preprocessing:
#   aligned with STAR, with given parameters
#
# DEPEND:
# bash(v4.2.46)
# gcc(v6.2.0)
# python(v3.6.0)
# java(v1.8.0_112)
# star(v2.5.4a)
# samtools(v1.3.1)
# picard(v2.8.0)
#
# Please, use --help|-h for help. 
# For more information, please, visit:
#       https://github.com/gimelbrantlab/ASE
#

import argparse
import os
import glob
import tempfile
import subprocess

def parse_param_lists(param, n_genomes):
    param_list = param.split(',')
    if (len(param_list)==1):
        return [param_list[0] for i in range(n_genomes)]
    elif (len(param_list)!=n_genomes):
        msg = "Incorrect list STAR_INDEX provided."
        raise argparse.ArgumentTypeError(msg)
    return param_list
 
def star_indexing(fa, gtf=None):
    if (gtf):
        cmd = ' '.join(["STAR --runMode genomeGenerate --genomeDir", os.path.dirname(fa), "--genomeFastaFiles", fa, "--sjdbGTFfile", gtf])
    else:
        cmd = ' '.join(["STAR --runMode genomeGenerate --genomeDir", os.path.dirname(fa), "--genomeFastaFiles", fa])
    print(cmd)
    subprocess.check_output(cmd, shell=True)
    return

def dict_to_string(D):
    return ' '.join([' '.join([item, D[item]]) for item in D])

def do_STAR(fastqs, genomes, odir, gtf=None, multimapNmax=None, just_names=False):
    results = []
    for genome in genomes:
        genome_idx = genomes.idx(genome)
        genome_name = os.path.splitext(os.path.basename(genome))[0]
        genome_mark = "_To" + genome_name + '.'

        cmd_sub = ' '.join(["--runThreadN 4 --outSAMtype BAM SortedByCoordinate", \
                           "--outSAMattrRGline", "ID:"+genome_name, \
                           "--genomeDir", os.path.dirname(genome)])
        if gtf:
            cmd_sub = ' '.join([cmd, "--sjdbGTFfile", gtf[genome_idx]])

        for fastq in fastqs:
            file_prefix = os.path.splitext(os.path.basename(fastq))[0] + genome_mark
            cmd = ' '.join(["STAR", "--readFilesIn", fastq, "--outFileNamePrefix", file_prefix, cmd_sub])
            if (multimapNmax):
                cmd = ' '.join([cmd, "--outFilterMultimapNmax", str(multimapNmax)])

            results += file_prefix + "Aligned.sortedByCoord.out.bam"
            if (just_names):
                continue
            print(cmd)
            subprocess.check_output(cmd, shell=True)
    return results

def do_PICARD_deduplication(bams):
    for bam in bams:
        # add SM tag in readgroups if needed:
        # samtools view -H $sample | sed -E 's/^@RG(.*)ID:(.*)/@RG\tID:\2\tSM:None/g' |  samtools reheader - $sample > $SAMPLEbam
        # Compose the read group identifier in the following format: @RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1     
        
        # coordinate sort of bam if needed:
        # java -jar $PICARD/picard-2.8.0.jar SortSam I=$SAMPLEbam O=$SAMPLESORTbam SORT_ORDER=coordinate

        imo = {'-I': bam, '-M': bam.replace(".bam", ".deduplicated.metrix"), '-O': bam.replace(".bam", ".deduplicated.bam")}
        cmd = ' '.join("gatk MarkDuplicates", dict_to_string(imo), "--REMOVE_DUPLICATES true")
        print(cmd)
        subprocess.check_output(cmd, shell=True)
    return

def main():
    # ARGUMENTS HANDLING: 
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", required=True, metavar="[method]", help="May be: CASTEL or ... .")
    parser.add_argument("--ALIGN", required=False, default="False", metavar="[True|False]", help="Boolean parameter for align need, eigther one for each genome or comma-separated list in the same order. If False will try to find corresponding bam in o_dir. Default: False.")
    parser.add_argument("--fastq_dir", required=True, metavar="[/path/to/dir]", help="Path to directory with fastq files (processing will be implemented fo all of them). Required.") 
    parser.add_argument("--o_dir", required=True, metavar="[/path/to/dir]", help="Path to output directory. Required.")
    parser.add_argument("--genome", required=True, metavar="[/path/to/fa1,/path/to/fa2]", help="Comma-separated list of fasta files, on which RNAseq will be alligned. Required.") 
    parser.add_argument("--gtf", required=False, metavar="[/path/to/file]", help="Path to reference gtf file(s). Should be one for each genome or comma-separated list in the same order.")
    parser.add_argument("--STAR_INDEX", required=False, default="False", metavar="[True|False]", help="Boolean parameter for star indexing need, eigther one for each genome or comma-separated list in the same order. Default: False.")
    parser.add_argument("--DEDUPLICATE", required=False, default="False", metavar="[True|False]", help="Boolean parameter for picard deduplication need, eigther one for each genome or comma-separated list in the same order. Default: False.")
    args = parser.parse_args()

    # SET VARS and TEST FOR CORRECT INPUT:
    fastq = glob.glob(os.path.join(args.fastq_dir, "*.fastq"))
    genome = args.STAR_INDEX.split(',')
    g_len = len(genome)
    if (args.gtf is not None):
        gtf = parse_param_lists(args.gtf, g_len)
    if (args.STAR_INDEX is not None):
        star_index = parse_param_lists(args.STAR_INDEX, g_len)
    if (args.DEDUPLICATE is not None):
        deduplicate = parse_param_lists(args.DEDUPLICATE, g_len)
    os.makedirs(args.o_dir, exist_ok=True)

    # IF STAR INDEX IS NEEDED:
    if (star_index):
        for i in range(g_len):
            if (star_index[i]=="True" and gtf):
                star_indexing(genome[i], gtf[i])
            elif (star_index[i]=="True" and not gtf):
                star_indexing(genome[i])

    # METHOD:
    if (args.method is not None):
        method = args.method
        if (method == "CASTEL"):
            bam = do_STAR(fastq, genome, args.o_dir, gtf, 1)
 
    # IF DEDUPLICATION:
    if (deduplicate):
        if (not bam):
            if (method =="CASTEL"):
                bam = do_STAR(fastq, genome, args.o_dir, gtf, 1, True)
        for i in range(g_len):
            if (deduplicte[i]=="True"):
                do_PICARD_deduplication(bam)

if __name__ == "__main__":
    main()

