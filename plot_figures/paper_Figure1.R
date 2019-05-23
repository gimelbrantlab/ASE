setwd("~/Dropbox (Partners HealthCare)/replicates_ASE/code/ASE/plot_figures")


library(tidyverse)
library(fitdistrplus)
library(cowplot)
library(ggrepel)

source("../R/ASE_functions.R")
source("../R/PerformAIAnalysis_CC.R")
source("../R/DownstreamASE.R")
source("../plot_figures/utilPlots.R")


# Load data
NEB_data_30mln = GetGatkPipelineTabs(paste0("../../../data/kidney/submln/", "MLN30_SG", 1:6, "_N955_", "NEB", "_R1_merged_v2.mln30_trial5_processed_gene_extended2.txt"), c(5,5,5,5,5,5), multiple = T)
SMART10ng_data_30mln = GetGatkPipelineTabs(paste0("../../../data/kidney/submln/", "MLN30_SG", 7:12, "_N955_", "SMARTseq_10_ng", "_R1_merged_v2.mln30_trial5_processed_gene_extended2.txt"), c(5,5,5,5,5,5), multiple = T)
SMART100pg_data_30mln = GetGatkPipelineTabs(paste0("../../../data/kidney/submln/", "MLN30_SG", 1:6, "_N955_", "SMARTseq_100_pg", "_R1_merged_v2.mln30_trial5_processed_gene_extended2.txt"), c(5,5,5,5,5,5), multiple = T)

data_30mln = list(NEB_data_30mln, SMART10ng_data_30mln, SMART100pg_data_30mln)
rm(NEB_data_30mln)
rm(SMART10ng_data_30mln)
rm(SMART100pg_data_30mln)

removeX <- function(DF, legitim_chrgenes){
  return(DF[DF$ensembl_gene_id %in% legitim_chrgenes$gene, ])
}
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

# Select pair of replicates for analysis (for seed=1 it is replicates 2 and 6)
set.seed(1)
sample2reps30mln = sample(0:5, 2, replace=F)*5 + sample(1:5, 2, replace=T)

sample2reps30mln_6 = sort((sample2reps30mln-1) %/% 5 + 1)
CCn_sample2reps30mln = sum((5:1)[0:(sample2reps30mln_6[1]-1)]) + sample2reps30mln_6[2]-sample2reps30mln_6[1]

CC_sample2_30mln = lapply(CC, function(x){x[CCn_sample2reps30mln]})

data_30mln_noX_sample2 = lapply(data_30mln_noX, function(x){
  x[, sort(c(1, sample2reps30mln*2, sample2reps30mln*2+1))]
})

# Get data frame with AIs and coverages for 2 selected replicates + tests outputs (T/F)

# Get data for NEB

DF_NEB_BT <- GetDataForExperiment_BT(df = data_30mln_noX_sample2, CC_df = CC_sample2_30mln, exp_name = "NEB", exp_n = 1)

# Plot figures with classification (middle panel)

figure_1C = ggplot(DF_NEB_BT[DF_NEB_BT$meanCOV.x < 9999, ], aes(meanCOV.x, AI.x)) +
  geom_point(aes(color=BT.x), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  facet_grid(~ BT.x) +
  theme_bw() +
  #ggtitle("Replicate 1") +
  labs(x = "Total Allelic Counts", y = "Gene AI - Replicate 1", color = "Binomial Test") +
  scale_x_continuous(trans='log10') +
  theme(legend.position="None",  text = element_text(size=18),
        strip.background = element_rect(fill="#E0E0E0"))

figure_1D = ggplot(DF_NEB_BT[DF_NEB_BT$meanCOV.y < 9999, ], aes(meanCOV.y, AI.y)) +
  geom_point(aes(color=BT.y), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  facet_grid(~ BT.x) +
  theme_bw() +
  #ggtitle("Replicate 2 relative to replicate 1") +
  labs(x = "Total Allelic Counts", y = "Gene AI - Replicate 2", color = "Binomial Test") +
  scale_x_continuous(trans='log10') +
  theme(legend.position="None", text = element_text(size=18),
        strip.background = element_rect(fill="#E0E0E0"))


# Repeat analysis for SMART-seq (right panel)

DF_SMART_10_BT <- GetDataForExperiment_BT(df = data_30mln_noX_sample2, CC_df = CC_sample2_30mln, exp_name = "SMARTseq 10ng", exp_n = 2)

figure_1E <- ggplot(DF_SMART_10_BT[DF_SMART_10_BT$meanCOV.y < 9999 & DF_SMART_10_BT$BT.x=="Imbalanced Genes \n (according to replicate 1)", ], aes(meanCOV.y, AI.y)) +
  geom_point(aes(color=BT.y), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  facet_grid(~ BT.x) +
  theme_bw() +
  labs(x = "Total Allelic Counts", y = "Gene AI  - Replicate 2", color = "Binomial Test") +
  scale_x_continuous(trans='log10') +
  theme(legend.position="None", text = element_text(size=18),
        strip.background = element_rect(fill="#E0E0E0"))

DF_SMART_0.1_BT <- GetDataForExperiment_BT(df = data_30mln_noX_sample2, CC_df = CC_sample2_30mln, exp_name = "SMARTseq 0.1ng", exp_n = 3)

figure_1F <- ggplot(DF_SMART_0.1_BT[DF_SMART_0.1_BT$meanCOV.y < 9999 & DF_SMART_0.1_BT$BT.x=="Imbalanced Genes \n (according to replicate 1)", ], aes(meanCOV.y, AI.y)) +
  geom_point(aes(color=BT.y), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  facet_grid(~ BT.x) +
  theme_bw() +
  labs(x = "Total Allelic Counts", y = "Gene AI  - Replicate 2", color = "Binomial Test") +
  scale_x_continuous(trans='log10') +
  theme(legend.position="None", text = element_text(size=18),
        strip.background = element_rect(fill="#E0E0E0"))

# Left panel

DF_NEB_BT$COV = rowMeans(DF_NEB_BT[, c("meanCOV.x", "meanCOV.y")])

figure_1A = ggplot(DF_NEB_BT, aes(AI.x, AI.y)) +
  geom_point(aes(color=COV), size=0.8) +
  scale_color_gradient(low="lightgray", high="black", trans = "log2",
                       breaks = c(10, 500, 5000),
                       labels = c(10, 500, 5000)) +
  theme_bw() +
  coord_fixed() +
  labs(x = "Gene AI - Replicate 1", y = "Gene AI - Replicate 2", color = "Allelic Counts") +
  theme(legend.key = element_blank(),legend.background=element_blank(),
        legend.position=c(0.2,0.8), text = element_text(size=18))

figure_1B = ggplot(DF_NEB_BT, aes(meanCOV.x, meanCOV.y)) +
  geom_point(color="black", size=0.8) +
  theme_bw() +
  coord_fixed() +
  labs(x = "Total Allelic Counts - Replicate 1", y = "Total Allelic Counts - Replicate 2") +
  theme(legend.position="bottom", text = element_text(size=18)) +
  xlim(0,10000) + ylim(0,10000)

PLT_fig1 = plot_grid(
  plot_grid(
    figure_1A,
    NULL,
    ncol=1,
    labels = c("A","C")
  ),
  plot_grid(
    figure_1B,
    NULL,
    ncol=1,
    labels = c("B","")
  ),
  plot_grid(
    figure_1C,
    figure_1D,
    ncol=1,
    labels = c("D","E")
  ),
  plot_grid(
    figure_1E,
    figure_1F,
    ncol=1,
    labels = c("F","G")
  ),
  ncol=4, scale = c(0.9, 0.9, 0.9, 0.9), rel_widths = c(1,1,1,0.55)
)

cowplot::save_plot(PLT_fig1, file="~/Dropbox (Partners HealthCare)/replicates_ASE/manuscript/Figures/fig.1_v2.pdf", base_height = 12, base_width = 24)


