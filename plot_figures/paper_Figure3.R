setwd("~/Dropbox (Partners HealthCare)/replicates_ASE/code/ASE/plot_figures")


library(tidyverse)
library(fitdistrplus)
library(cowplot)
library(ggrepel)
library(wesanderson)

source("../R/ASE_functions.R")
source("../R/PerformAIAnalysis_CC.R")
source("../R/DownstreamASE.R")
source("../plot_figures/utilPlots.R")

methods = c("Experiment 1", "Experiment 2", "Experiment 3")


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

Ntrials = 1
set.seed(1)
sampleMreps = sapply(1:(max(Ntrials, 10)), function(i){(sample(1:6, 2)-1)*5 + sample(5, 2, replace=T)})
list_of_datas = data_30mln_noX
list_of_consts = list(mean(CC[[1]]), mean(CC[[2]]), mean(CC[[3]]))
list_of_libprepnames = list("1: NEBNext (100ng)", "2: SMARTseq (10ng)", "3: SMARTseq (0.1ng)")
reppair = sampleMreps[, 1]


DF_forplot = CreateForplotDF_btNbtcc(lapply(1:3, function(i){
  CreateForplotDF(list_of_datas[[i]], reppair, list_of_consts[[i]], list_of_libprepnames[[i]])
}))

DF_forplot$BTBFCC$BT.x = as.factor(DF_forplot$BTBFCC$BT.x)
levels(DF_forplot$BTBFCC$BT.x) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")
DF_forplot$BTBFCC$BT.y = as.factor(DF_forplot$BTBFCC$BT.y)
levels(DF_forplot$BTBFCC$BT.y) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")


# Calculate disagreement percentages

# For before correction

pdis = lapply(1:3, function(m){
  sapply(1:ncol(combn(6, 2)), function(j){
    sample2reps30mln =(combn(6, 2)[, j]-1)*5 + sample(1:5, 2, replace=T)

    sample2reps30mln_6 = sort((sample2reps30mln-1) %/% 5 + 1)
    CCn_sample2reps30mln = sum((5:1)[0:(sample2reps30mln_6[1]-1)]) + sample2reps30mln_6[2]-sample2reps30mln_6[1]

    CC_sample2_30mln = CC[[m]][CCn_sample2reps30mln]

    data_30mln_sample2 = data_30mln[[m]][, sort(c(1, sample2reps30mln*2, sample2reps30mln*2+1))]

    p = GetPercDiffForExperiment_BT(df=data_30mln_sample2, CC_df=CC_sample2_30mln, exp_name = "any")

    return(p)
  })
})

# And after correction

pdis_CC = lapply(1:3, function(m){
  sapply(1:ncol(combn(6, 2)), function(j){
    sample2reps30mln =(combn(6, 2)[, j]-1)*5 + sample(1:5, 2, replace=T)

    sample2reps30mln_6 = sort((sample2reps30mln-1) %/% 5 + 1)
    CCn_sample2reps30mln = sum((5:1)[0:(sample2reps30mln_6[1]-1)]) + sample2reps30mln_6[2]-sample2reps30mln_6[1]

    CC_sample2_30mln = CC[[m]][CCn_sample2reps30mln]

    data_30mln_sample2 = data_30mln[[m]][, sort(c(1, sample2reps30mln*2, sample2reps30mln*2+1))]

    p = GetPercDiffForExperiment_BT_CC(df=data_30mln_sample2, CC_df=CC_sample2_30mln, exp_name="any")

    return(p)
  })
})

df = rbind(
  do.call(rbind, lapply(1:3, function(i){data.frame(P=as.vector(pdis[[i]]), libprep=methods[i], what="before correction")})),
  do.call(rbind, lapply(1:3, function(i){data.frame(P=as.vector(pdis_CC[[i]]), libprep=methods[i], what="corrected")}))
)

# Figure 3C

res62_df <-
  lapply(1:3, function(m){
    df <- sapply(1:ncol(combn(6, 2)), function(j){
      pairs_sample = (combn(6, 2)[, j]-1)*5 + sample(1:5, 2, replace=T)
      for_6_df <- CountsToAI(data_30mln[[m]], meth="mergedToProportion")$AI
      res62 <- PerformBinTestAIAnalysisForConditionNPointVect_knownCC(data_30mln[[m]], vectReps=pairs_sample,
                                                                      vectRepsCombsCC = CC[[m]],
                                                                      ptVect = for_6_df,
                                                                      thr=10)
      FP_BTCC <- sum(res62$BT_CC, na.rm = T)
      FP_BT <- sum(res62$BT, na.rm = T)
      all_BTCC <- length(res62$BT_CC) - sum(is.na(res62$BT_CC))
      all_BT <- length(res62$BT) - sum(is.na(res62$BT))
      return(c(FP_BTCC, FP_BT, all_BTCC, all_BT))
    })
    df <- data.frame(t(df))
    colnames(df) <- c("FP_BTCC", "FP_BT", "all_BTCC", "all_BT")
    df$FP_BTCC_rate <- df$FP_BTCC / df$all_BTCC
    df$FP_BT_rate <- df$FP_BT / df$all_BT
    return(df)
  })


# res62_df <- data.frame(t(res62_df))
# colnames(res62_df) <- c("FP_BTCC", "FP_BT", "all_BTCC", "all_BT")
# res62_df$FP_BTCC_rate <- res62_df$FP_BTCC / res62_df$all_BTCC
# res62_df$FP_BT_rate <- res62_df$FP_BT / res62_df$all_BT

res62_df_NEB <- res62_df[[1]]
res62_df_NEB$experiment <- "NEB"

res62_df_SMART10 <- res62_df[[2]]
res62_df_SMART10$experiment <- "SMART10"

res62_df_SMART100 <- res62_df[[3]]
res62_df_SMART100$experiment <- "SMART100"

res62_df_all <- data.frame(rbind(res62_df_NEB, res62_df_SMART10, res62_df_SMART100))
res62_df_all <- res62_df_all[,c(7,5,6)] %>% gather(key="method", value = "FP_rate", -experiment)

# Plotting

figure_3A <- ggplot(DF_forplot$BTBFCC[DF_forplot$BTBFCC$meanCOV.y<9999,], aes(meanCOV.y, AI.y)) +
  geom_point(aes(color=BT.y), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  facet_grid(libprep ~ BT.x) +
  theme_bw() +
  labs(x = "Total Allele Counts", y = "Gene AI - Replicate 2", color = "Binomial Test") +
  scale_x_continuous(trans='log2') +
  theme(legend.position="None", text = element_text(size=18))

figure_3B <- ggplot(df, aes(x=libprep, y=P, col=libprep)) +
  geom_boxplot() +
  facet_grid(~what) +
  xlab("") +
  ylab("Disagreement % \nbetween two replicates") +
  theme_bw() +
  theme(legend.position="bottom", text = element_text(size=18)) +
  guides(col=guide_legend(title="", nrow = 1)) +
  scale_color_manual(labels=methods,
                     values=wes_palette(n=3, name="FantasticFox1")) +
  #scale_color_manual(values=c("#999999", "#66FFB2"), labels=c("no correction", "QCC correction")) +
  scale_x_discrete(labels=c("", "", ""))

figure_3D <- ggplot(res62_df_all, aes(x=method, y=FP_rate, col=experiment)) +
  geom_boxplot() +
  xlab("") +
  ylab("False positive rate") +
  theme_bw() +
  theme(legend.position="None", text = element_text(size=18)) +
  scale_color_manual(labels=methods,
                     values=wes_palette(n=3, name="FantasticFox1")) +
  scale_x_discrete(labels=c("no correction", "QCC correction"))

PLT_fig3 = plot_grid(
  plot_grid(figure_3A, figure_3B, labels = c("A", "B"), rel_widths = c(0.8, 0.8)),
  plot_grid(NULL, figure_3D, labels = c("C", "D"), rel_widths = c(0.8, 0.8)),
  nrow = 2,
  scale = c(0.9, 0.9, 0.9, 0.9)
)

cowplot::save_plot(PLT_fig3, file="~/Dropbox (Partners HealthCare)/replicates_ASE/manuscript/Figures/fig.3_v0.pdf", base_height = 12, base_width = 14)
