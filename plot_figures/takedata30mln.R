library(tidyverse)
library(fitdistrplus)
library(cowplot)
library(ggrepel)
library(wesanderson)

methods = c("Experiment 1", "Experiment 2", "Experiment 3")

# Load data
NEB_data_30mln = GetGatkPipelineTabs(paste0("../../../data/kidney/submln/", "MLN30_SG", 1:6, "_N955_", "NEB", "_R1_merged_v2.mln30_trial5_processed_gene_extended2.txt"), c(5,5,5,5,5,5), multiple = T)
SMART10ng_data_30mln = GetGatkPipelineTabs(paste0("../../../data/kidney/submln/", "MLN30_SG", 7:12, "_N955_", "SMARTseq_10_ng", "_R1_merged_v2.mln30_trial5_processed_gene_extended2.txt"), c(5,5,5,5,5,5), multiple = T)
SMART100pg_data_30mln = GetGatkPipelineTabs(paste0("../../../data/kidney/submln/", "MLN30_SG", 1:6, "_N955_", "SMARTseq_100_pg", "_R1_merged_v2.mln30_trial5_processed_gene_extended2.txt"), c(5,5,5,5,5,5), multiple = T)

data_30mln = list(NEB_data_30mln, SMART10ng_data_30mln, SMART100pg_data_30mln)
rm(NEB_data_30mln)
rm(SMART10ng_data_30mln)
rm(SMART100pg_data_30mln)

chrgenes <- read.delim("../../../data/kidney/submln/MLN30_SG3_N955_NEB_R1_merged_v2.mln30_trial5_processed_gene_extended2.txt")[, c("ensembl_gene_id", "chr")]
data_30mln_noX = lapply(data_30mln, function(x){
  x[x$ensembl_gene_id %in% chrgenes[chrgenes$chr!="chrX" & chrgenes$chr!="chrY", "ensembl_gene_id"], ]
})

# Load pre-calculated CCs
CC = lapply(c("../../../data/kidney/submln/NEB_CorrConsts_30mln_1.05.RData",
              "../../../data/kidney/submln/SMARTseq_10_ng_CorrConsts_30mln_1.05.RData",
              "../../../data/kidney/submln/SMARTseq_100_pg_CorrConsts_30mln_1.05.RData"),
            function(file){
              load(file)
              sapply(out_XXmln_SMART10ng, function(x){x$fittedCC})
            })
