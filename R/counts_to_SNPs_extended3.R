#!/usr/bin/env Rscript
# ***********************************************
# Title       : Allelic Imbalance Pipeline
# Description : Script that takes $name_stat_final.txt files as input,
#               merges these files, intersects counts with exons annotaion
#               
# Author      : Svetlana Vinogradova
# Date        : 08/24/17
# ***********************************************
#
# run:
#
# Rscript --vanilla counts_to_SNPs_extended2.R -p "/Users/svetlana/Dropbox (Partners HealthCare)/ASE/resultsUpdAugust/Abelson" -r Aleson_test -n Abl1,Abl2 -s GRCm38 -e /Users/svetlana -t exon -m /Users/svetlana/extract_vcf.txt
#
#
library("optparse", lib.loc="/home/sv111/R-3.5.1/library")
library("crayon", lib.loc="/home/sv111/R-3.5.1/library")
library("tidyverse", lib.loc="/home/sv111/R-3.5.1/library")
#
option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path to directory contatining input_data folder with allele counts", metavar="character"),
  make_option(c("-r", "--pr_name"), type="character", default="clones", 
              help="project name, for example Abelson_clones [default= %default]", metavar="character"),
  make_option(c("-s", "--sp_name"), type="character", default="GRCm38", 
              help="organism name, GRCm38 (mm10) or GRCh38 (hg38) [default= %default]", metavar="character"),
  make_option(c("-n", "--names"), type="character", default=NULL, 
              help="names of clones in comma-separated format, for example clone1,clone2,clone3 ", metavar="character"),
  make_option(c("-e", "--exons_path"), type="character", default=NULL, 
              help="path to directory with *exons.bed and *exons_full.txt files", metavar="character"),
  make_option(c("-t", "--type"), type="character", default="exon",
              help="exon or intron -- for exonic and intronic SNPs", metavar="character"),
  make_option(c("-m", "--snp_tab"), type="character", default=NULL,
              help="path to table with SNPs (from vcf): positions, names, info", metavar="character"),
  make_option(c("-a", "--actual_name"), type="character", default=".stat_0.txt",
	      help="provide the name suffix for the input files", metavar="character")
); 
#
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
# Check that all arguments are passed
if ((is.null(opt$path))|(is.null(opt$names))|(is.null(opt$exons_path))|(is.null(opt$snp_tab))){
  print_help(opt_parser)
  stop("Paths and names must be supplied", call.=FALSE)
}
# Find paths to exons
if (opt$sp_name %in% c("GRCm38", "mm10")) {
  exon_full <- paste0(opt$exons_path, "/GRCm38.v96.exons_full.txt")
  exon_bed <-  paste0(opt$exons_path, "/GRCm38.v96.exons.bed")
} else if (opt$sp_name %in% c("GRCh38", "hg38")) {
  exon_full <-  paste0(opt$exons_path, "/GRCm38.v96.exons_full.txt")
  exon_bed <-  paste0(opt$exons_path, "/GRCm38.v96.exons.bed")
} else {
  print_help(opt_parser)
  stop("Unknown organism", call.=FALSE)
}
#
# Read files created by python allele_counter
setwd(opt$path)
names_list <- unlist(strsplit(opt$names, ","))
# Read data for SNP positions
merged_counts <- read_delim(opt$snp_tab, delim="\t", escape_double = FALSE, trim_ws = TRUE, col_types = cols(`#CHROM` = col_character()))
colnames(merged_counts) <- c("contig" ,"position", "variantID", "refAllele", "altAllele")

print(unique(merged_counts$contig))

# Read data for all clones and merge together
for (count in (1:length(names_list))) {
  name <- paste0("input_data/", names_list[count], opt$actual_name)
  rep_tab <- read_table2(name, col_types = cols(contig = col_character()))  
  #rep_tab <- read_delim(name, " ", escape_double = FALSE, col_types = cols(contig = col_character()), trim_ws = TRUE)
  rep_tab <- rep_tab[ ,c("contig" ,"position", "refCount", "altCount")]
  merged_counts <- merge(merged_counts, rep_tab, by.x=c("contig", "position"), by.y=c("contig", "position"), all.x = TRUE)
  colnames(merged_counts)[1:3] <- c("contig" ,"position", "variantID")
}
merged_counts[is.na(merged_counts)] <- 0
al_names <- c("contig" ,"position", "variantID", "refAllele", "altAllele")
names_list_write <- paste0("rep", 1:length(names_list))
for (count in (1:length(names_list))) {
  al_names <- c(al_names, paste("ref", names_list_write[count], sep = "_"), paste("alt", names_list_write[count], sep = "_")) 
}
colnames(merged_counts) <- al_names
#
print(unique(merged_counts$contig))
# Create bed file for bedtools
merged_counts$chr <- paste0("chr",merged_counts$contig)
merged_counts$unique_id <- paste0(merged_counts$contig, merged_counts$position)
merged_counts$end <- merged_counts$position+1
merged_counts_bed <- merged_counts[,c('chr','position', 'end', 'unique_id')]
bed_name <- paste0(opt$pr_name, "_merged_extended2.bed")
bed_name_exons <- paste0(opt$pr_name, "_merged_to_exons_extended2.bed")
write.table(merged_counts_bed, file=bed_name, quote = F, row.names = F, col.names = F, sep="\t")
# Run bedtools
cmd <- paste0("/home/sv111/tools/bedtools2/bin/bedtools intersect -b ", exon_bed)
cmd <- paste0(cmd, " -a ")
cmd <- paste0(cmd, bed_name)
cmd <- paste0(cmd, " -wa -wb > ")
cmd <- paste0(cmd, bed_name_exons)
system(cmd)
# Read bedtools output
to_exons <- read_delim(bed_name_exons, delim="\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE, col_types = "ciicciic")
# Read exons full information
exons <- read_delim(exon_full, delim="\t", escape_double = FALSE, trim_ws = TRUE, col_names = T, col_types = cols(chromosome_name = col_character()))
#
# Merge back to exon annotation and full counts information
to_exons <- merge(to_exons, exons, by.x="X8", by.y="unique_id")
to_exons <- merge(to_exons, merged_counts, by.x="X4", by.y="unique_id")
to_exons <- to_exons[,c(22:(29+((length(names_list)-1)*2)), 9,10, 12:21)]
colnames(to_exons)[1] <- "chr"
to_exons <- to_exons %>% 
  select(-2)
merged_name <- paste0(opt$pr_name, "_merged_counts_extended2.txt")
write.table(to_exons, file=merged_name, quote = F, row.names = F, col.names = T, sep="\t")

