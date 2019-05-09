
library(tidyverse)
library(fitdistrplus)
library(cowplot)
library(ggrepel)

source("../R/ASE_functions.R")
source("../R/PerformAIAnalysis_CC.R")
source("../R/DownstreamASE.R")


# Load data
SMART10ng_data_20mln = GetGatkPipelineTabs(paste0("../../../data/kidney/submln/", "MLN20_SG", 7:12, "_N955_", "SMARTseq_10_ng", "_R1_merged_v2.mln20_trial5_processed_gene_extended2.txt"), c(5,5,5,5,5,5), multiple = T)
SMART100pg_data_20mln = GetGatkPipelineTabs(paste0("../../../data/kidney/submln/", "MLN20_SG", 1:6, "_N955_", "SMARTseq_100_pg", "_R1_merged_v2.mln20_trial5_processed_gene_extended2.txt"), c(5,5,5,5,5,5), multiple = T)
NEB_data_20mln = GetGatkPipelineTabs(paste0("../../../data/kidney/submln/", "MLN20_SG", 1:6, "_N955_", "NEB", "_R1_merged_v2.mln20_trial5_processed_gene_extended2.txt"), c(5,5,5,5,5,5), multiple = T)

data_20mln = list(NEB_data_20mln, SMART10ng_data_20mln, SMART100pg_data_20mln)

removeX <- function(DF, legitim_chrgenes){
  return(DF[DF$ensembl_gene_id %in% legitim_chrgenes$gene, ])
}
chrgenes <- read.delim("../../../data/kidney/submln/MLN20_SG3_N955_NEB_R1_merged_v2.mln20_trial5_processed_gene_extended2.txt")[, c("ensembl_gene_id", "chr")]

data_20mln_noX = lapply(data_20mln, function(x){
  x[x$ensembl_gene_id %in% chrgenes[chrgenes$chr!="chrX" & chrgenes$chr!="chrY", "ensembl_gene_id"], ]
})


# Load pre-calculated CCs
CC = list(
  c(2.342191,2.265159,2.402253,2.352722,2.375676,1.810909,1.984972,1.894131,1.936703,1.919385,1.850982,1.892572,1.923280,1.975188,1.886451),
  c(1.591385,1.621971,1.640987,1.587173,1.617223,1.571266,1.610344,1.555469,1.631672,1.580222,1.611677,1.596801,1.627877,1.605522,1.645494),       
  c(2.905560,2.818407,2.830344,2.754865,2.755437,2.912273,2.949068,2.858991,2.896754,2.834602,2.844827,2.892253,2.803532,2.869579,2.852329)
)

# Select pair of replicates for analysis (for seed=1 it is replicates 2 and 6)
set.seed(1)
sample2reps20mln = sample(0:5, 2, replace=F)*5 + sample(1:5, 2, replace=T)

sample2reps20mln_6 = sort((sample2reps20mln-1) %/% 5 + 1)
CCn_sample2reps20mln = sum((5:1)[0:(sample2reps20mln_6[1]-1)]) + sample2reps20mln_6[2]-sample2reps20mln_6[1]

CC_sample2_20mln = lapply(CC, function(x){x[CCn_sample2reps20mln]})

data_20mln_noX_sample2 = lapply(data_20mln_noX, function(x){
  x[, sort(c(1, sample2reps20mln*2, sample2reps20mln*2+1))]
})

# Get data frame with AIs and coverages for 2 selected replicates + tests outputs (T/F)
forbtplots1A = CreateForplotDF(data_20mln_noX_sample2[[1]], 1:2, 
                               pairconst=CC_sample2_20mln[[1]], 
                               "NEB")

# Split two tests outputs in separate datasets 
DF_forplot_1A = CreateForplotDF_btNbtcc(list(forbtplots1A))


# Calculate % of disagreement between replicates
percent_of_diff_color_df_1A = CreateForplotDF_btNbtcc_colorescapers(DF_forplot_1A)
percent_of_diff_color_df_1A$BTBF$test = "binomial"
percent_of_diff_color_df_1A_BT = percent_of_diff_color_df_1A$BTBF

percent_of_diff_color_df_1A_BT$BT.x = as.factor(percent_of_diff_color_df_1A_BT$BT.x)
levels(percent_of_diff_color_df_1A_BT$BT.x) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")

# Get data for Binomial test + Bonferroni correction only
DF_forplot_1A_BT = DF_forplot_1A$BTBF
DF_forplot_1A_BT$BT.x = as.factor(DF_forplot_1A_BT$BT.x)
levels(DF_forplot_1A_BT$BT.x) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)") 
DF_forplot_1A_BT$BT.y = as.factor(DF_forplot_1A_BT$BT.y)
levels(DF_forplot_1A_BT$BT.y) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)") 


# Plot figures with classification 
figure_1C = ggplot(DF_forplot_1A_BT[DF_forplot_1A_BT$meanCOV.x < 9999, ], aes(meanCOV.x, AI.x)) +
  geom_point(aes(color=BT.x), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  facet_grid(~ BT.x) +
  theme_bw() +
  #ggtitle("Replicate 1") +
  labs(x = "Total Allelic Counts", y = "Gene AI - Replicate 1", color = "Binomial Test") +
  scale_x_continuous(trans='log10') +
  theme(legend.position="None",  text = element_text(size=18)) 

figure_1D = ggplot(DF_forplot_1A_BT[DF_forplot_1A_BT$meanCOV.y < 9999, ], aes(meanCOV.y, AI.y)) +
  geom_point(aes(color=BT.y), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  facet_grid(~ BT.x) +
  theme_bw() +
  #ggtitle("Replicate 2 relative to replicate 1") +
  labs(x = "Total Allelic Counts", y = "Gene AI - Replicate 2", color = "Binomial Test") +
  scale_x_continuous(trans='log10') +
  theme(legend.position="None", text = element_text(size=18)) 


# Repeat analysis for SMART-seq

forbtplotsDF1C = lapply(1:3, function(i){
  forbtplots1C = CreateForplotDF(data_20mln_noX_sample2[[i]], 1:2, 
                                 pairconst=CC_sample2_20mln[[i]], 
                                 "SMARTseq 0.1ng")
  DF_forplot_1C = CreateForplotDF_btNbtcc(list(forbtplots1C))
  
  percent_of_diff_color_df_1C = CreateForplotDF_btNbtcc_colorescapers(DF_forplot_1C)
  percent_of_diff_color_df_1C$BTBF$test = "binomial"
  percent_of_diff_color_df_1C_BT = percent_of_diff_color_df_1C$BTBF
  
  percent_of_diff_color_df_1C_BT$BT.x = as.factor(percent_of_diff_color_df_1C_BT$BT.x)
  levels(percent_of_diff_color_df_1C_BT$BT.x) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")
  
  
  DF_forplot_1C_BT = DF_forplot_1C$BTBF
  DF_forplot_1C_BT$BT.x = as.factor(DF_forplot_1C_BT$BT.x)
  levels(DF_forplot_1C_BT$BT.x) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)") 
  DF_forplot_1C_BT$BT.y = as.factor(DF_forplot_1C_BT$BT.y)
  levels(DF_forplot_1C_BT$BT.y) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)") 
  
  return(list(percent_of_diff_color_df_1C_BT, DF_forplot_1C_BT, DF_forplot_1C$BTBF))  
})

figure_1E = lapply(1:2, function(i){
  
  df = forbtplotsDF1C[[i]][[2]]
  pcnt = forbtplotsDF1C[[i]][[1]]
  
  plt_bt_1C = ggplot(df[df$meanCOV.y < 9999 & df$BT.x=="Imbalanced Genes \n (according to replicate 1)", ], aes(meanCOV.y, AI.y)) +
    geom_point(aes(color=BT.y), size=0.8) +
    scale_color_manual(values=c("salmon","royalblue1")) +
    facet_grid(~ BT.x) +
    theme_bw() +
    labs(x = "Total Allelic Counts - Replicate 2", y = "Gene AI  - Replicate 2", color = "Binomial Test") +
    scale_x_continuous(trans='log10') +
    theme(legend.position="None", text = element_text(size=18)) 
  return(plt_bt_1C)
})

figure_1A = lapply(1:2, function(i){
  
  df = forbtplotsDF1C[[i]][[3]]
  df$COV = rowMeans(df[, c("meanCOV.x", "meanCOV.y")])
  
  plt_bt_1C = ggplot(df, aes(AI.x, AI.y)) +
    geom_point(aes(color=COV), size=0.8) +
    scale_color_gradient(low="lightgray", high="black", trans = "log2", 
                         breaks = c(10, 500, 5000),
                         labels = c(10, 500, 5000)) +
    theme_bw() +
    coord_fixed() +
    labs(x = "Gene AI - Replicate 1", y = "Gene AI - Replicate 2", color = "Total Allele Counts") +
    theme(legend.key = element_blank(),legend.background=element_blank(),
          legend.position=c(0.3,0.8), text = element_text(size=18)) 
  return(plt_bt_1C)
})

figure_1B = lapply(1:3, function(i){
  
  df = forbtplotsDF1C[[i]][[3]]
  
  plt_bt_1C = ggplot(df, aes(meanCOV.x, meanCOV.y)) +
    geom_point(color="black", size=0.8) +
    theme_bw() +
    coord_fixed() +
    labs(x = "Total Allelic Counts - Replicate 1", y = "Total Allelic Counts - Replicate 2") +
    theme(legend.position="bottom", text = element_text(size=18)) +
    xlim(0,10000) + ylim(0,10000)
  
  return(plt_bt_1C)
})

PLT_fig1 = plot_grid(
  plot_grid(
    figure_1A[[1]],
    figure_1B[[1]],
    ncol=1, 
    labels = c("A","B")
  ),
  plot_grid(
    figure_1C,
    figure_1D,
    ncol=1, 
    labels = c("C","D")
  ),
  plot_grid(
    figure_1E[[1]],
    figure_1E[[2]],
    ncol=1, 
    labels = c("E","F")
  ),
  ncol=3, scale = c(0.9, 0.9, 0.9), rel_widths = c(1,1,0.55)
)


cowplot::save_plot(PLT_fig1, file="~/Dropbox (Partners HealthCare)/replicates_ASE/manuscript/Figures/fig.1_v1.2.pdf", base_height = 12, base_width = 18)


