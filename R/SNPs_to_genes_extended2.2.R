# ***********************************************
# Title       : Allelic Imbalance Pipeline
# Description : Script that takes $name_merged_counts.txt files as input,
#               merges information from SNPs, 
#               outputs $name_processed.txt
#               
# Author      : Svetlana Vinogradova (+ Asya Mendelevich)
# Date        : 03/28/18
# ***********************************************
#
#
#
# run:
#
# Rscript SNPs_to_genes.R -p "/n/scratch2/sv111/data/kidney/results/input_data/NEB_100ng/" -n NEB_100ng/ -k 6
#
#
# Libraries needed
library("optparse", lib.loc="/home/sv111/R-3.5.1/library")
library("crayon", lib.loc="/home/sv111/R-3.5.1/library")
library("tidyverse", lib.loc="/home/sv111/R-3.5.1/library")
#
option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="path to directory contatining merged_counts file", metavar="character"),
  make_option(c("-n", "--file_name"), type="character", default=NULL, 
              help="name of the experiment", metavar="character"),
  make_option(c("-k", "--n"), type="integer", default=2, 
              help="number of samples/replicates [default= %default]")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
n <- opt$n
file_in <- paste0(opt$path, paste0(opt$file_name, "_merged_counts_extended2.txt"))
file_out_g <- paste0(opt$path, paste0(opt$file_name, "_processed_gene_extended2.txt"))
file_out_s <- paste0(opt$path, paste0(opt$file_name, "_processed_snp_extended2.txt"))

# Read file and deduplicate rows for each pair SNP-GENE:
df <- read_delim(file_in, "\t", escape_double = FALSE, trim_ws = TRUE)
df <- df[,c(1:(6+2*n))]
df <- add_count(df, df$chr, df$position)
df <- df[!duplicated(df),]

# Prepare SNP table:
snp_df <- data.frame(snp_id = paste(df$chr, df$position, df$ensembl_gene_id,df$variantID, sep = '_'),
                     df[,c(6:(5+(n*2)))], aggr_n = df$n, chr = df$chr, pos = df$position, gene_id = df$ensembl_gene_id)
snp_df <- snp_df[order(snp_df$chr, snp_df$pos), ]

# Aggregate by gene name: sum and concatenate:
names <- paste0("rep", 1:n)
aggr <- aggregate(df[,c(6:(5+(n*2)))], by=list(df$ensembl_gene_id), sum)
colnames(aggr)[1] <- "ensembl_gene_id"

aggr_coll <- aggregate(df[,c(6:(5+(n*2)),dim(df)[2],1,2)], by=list(df$ensembl_gene_id), paste, collapse = ", ")
agg_names_ref <- "ensembl_gene_id"
for (count in (1:n)) {
  agg_names_ref <- c(agg_names_ref, paste("aggr_ref", names[count], sep = "_"),paste("aggr_alt", names[count], sep = "_")) 
}
agg_names_ref <- c(agg_names_ref, "aggr_n", "chr", "pos")
colnames(aggr_coll) <- agg_names_ref
aggr_coll$chr <- sapply(aggr_coll$chr, function(x){strsplit(x,',')[[1]][1]})

# Merge all together:
aggr_merged <- unique(merge(aggr, aggr_coll, by.x = "ensembl_gene_id", by.y="ensembl_gene_id"))
aggr_merged <- aggr_merged[order(aggr_merged$chr, sapply(aggr_merged$pos, function(x){as.integer(strsplit(x,',')[[1]][1])})), ]

# Write the final tables:
write.table(aggr_merged, file_out_g, sep="\t", quote = F, row.names = F)
write.table(snp_df, file_out_s, sep="\t", quote = F, row.names = F)

