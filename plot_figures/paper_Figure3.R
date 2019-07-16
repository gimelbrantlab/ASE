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

methods = c("NEBNext (100ng)", "SMARTseq (10ng)", "SMARTseq (0.1ng)")


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
list_of_libprepnames = list("NEBNext (100ng)", "SMARTseq (10ng)", "SMARTseq (0.1ng)")
reppair = sampleMreps[, 1]

CC_df <- data.frame(t(rbind(CC[[1]], CC[[2]], CC[[3]])))
colnames(CC_df) <- c("NEBNext (100ng)", "SMARTseq (10ng)", "SMARTseq (0.1ng)")
CC_DF <- reshape2::melt(CC_df)


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

df$libprep <- as.factor(df$libprep)
df$libprep <- factor(df$libprep, levels=c("NEBNext (100ng)","SMARTseq (10ng)","SMARTseq (0.1ng)"))


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



###########

# Circles

# Functions

CalculateTestWithHalf1 <- function(df, cc, exp_name, thr=10){
  res = Reduce(function(x, y) merge(x, y, by="ID"),
               lapply(1:(ncol(df)%/%(2)), function(j){
                 PerformBinTestAIAnalysisForConditionNPoint_knownCC(df, j, cc, thr=thr)[, c("ID", "BT", "BT_CC")]
               })
  )
  names(res) = c("ID",
                 paste0(c("BT.", "BT_CC."),
                        rep(1:(ncol(df)%/%(2)), each=2)
                 )
  )
  return(res)
}

# Data prep

#NEB

## single rep
df_NEB_2_3 = CalculateTestWithHalf1(data_30mln_noX[[1]][, c(1, 2*5*2, 2*5*2+1, 3*5*2, 3*5*2+1)],
                                    CC[[1]][6], list_of_libprepnames[[1]], thr=10)
## merge two reps
df_NEB_24 = CalculateTestWithHalf1(MergeSumCounts(data_30mln_noX[[1]][, c(1, 2*5*2, 2*5*2+1, 4*5*2, 4*5*2+1)]),
                                   CC[[1]][7], list_of_libprepnames[[1]], thr=10)
df_NEB_35 = CalculateTestWithHalf1(MergeSumCounts(data_30mln_noX[[1]][, c(1, 3*5*2, 3*5*2+1, 5*5*2, 5*5*2+1)]),
                                   CC[[1]][11], list_of_libprepnames[[1]], thr=10)
df_NEB_24_35 = merge(df_NEB_24, df_NEB_35, by="ID")


#SMART100

## single rep
df_SM100_2_3 = CalculateTestWithHalf1(data_30mln_noX[[3]][, c(1, 2*5*2, 2*5*2+1, 3*5*2, 3*5*2+1)],
                                      CC[[3]][6], list_of_libprepnames[[3]], thr=10)
## merge two reps
df_SM100_24 = CalculateTestWithHalf1(MergeSumCounts(data_30mln_noX[[3]][, c(1, 2*5*2, 2*5*2+1, 4*5*2, 4*5*2+1)]),
                                     CC[[3]][7], list_of_libprepnames[[3]], thr=10)
df_SM100_35 = CalculateTestWithHalf1(MergeSumCounts(data_30mln_noX[[3]][, c(1, 3*5*2, 3*5*2+1, 5*5*2, 5*5*2+1)]),
                                     CC[[3]][11], list_of_libprepnames[[3]], thr=10)
df_SM100_24_35 = merge(df_SM100_24, df_SM100_35, by="ID")


#SMART10

## single rep
df_SM10_2_3 = CalculateTestWithHalf1(data_30mln_noX[[2]][, c(1, 2*5*2, 2*5*2+1, 3*5*2, 3*5*2+1)],
                                     CC[[2]][6], list_of_libprepnames[[2]], thr=10)
## merge two reps
df_SM10_24 = CalculateTestWithHalf1(MergeSumCounts(data_30mln_noX[[2]][, c(1, 2*5*2, 2*5*2+1, 4*5*2, 4*5*2+1)]),
                                    CC[[2]][7], list_of_libprepnames[[2]], thr=10)
df_SM10_35 = CalculateTestWithHalf1(MergeSumCounts(data_30mln_noX[[2]][, c(1, 3*5*2, 3*5*2+1, 5*5*2, 5*5*2+1)]),
                                    CC[[2]][11], list_of_libprepnames[[2]], thr=10)
df_SM10_24_35 = merge(df_SM10_24, df_SM10_35, by="ID")

# Euler

eul_NEB_2_3 = euler(na.omit(df_NEB_2_3[, c("BT.1","BT.2")]), shape='ellipse')
eul_NEB_24_35_bt = euler(na.omit(df_NEB_24_35[, c("BT.1.x","BT.1.y")]), shape='ellipse')
eul_NEB_24_35_btcc = euler(na.omit(df_NEB_24_35[, c("BT_CC.1.x","BT_CC.1.y")]), shape='ellipse')

NEB_eul_list <- list(eul_NEB_2_3, eul_NEB_24_35_bt, eul_NEB_24_35_btcc)

eul_SM100_2_3 = euler(na.omit(df_SM100_2_3[, c("BT.1","BT.2")]), shape='ellipse')
eul_SM100_24_35_bt = euler(na.omit(df_SM100_24_35[, c("BT.1.x","BT.1.y")]), shape='ellipse')
eul_SM100_24_35_btcc = euler(na.omit(df_SM100_24_35[, c("BT_CC.1.x","BT_CC.1.y")]), shape='ellipse')

SMART100_eul_list <- list(eul_SM100_2_3, eul_SM100_24_35_bt, eul_SM100_24_35_btcc)

eul_SM10_2_3 = euler(na.omit(df_SM10_2_3[, c("BT.1","BT.2")]), shape='ellipse')
eul_SM10_24_35_bt = euler(na.omit(df_SM10_24_35[, c("BT.1.x","BT.1.y")]), shape='ellipse')
eul_SM10_24_35_btcc = euler(na.omit(df_SM10_24_35[, c("BT_CC.1.x","BT_CC.1.y")]), shape='ellipse')

SMART10_eul_list <- list(eul_SM10_2_3, eul_SM10_24_35_bt, eul_SM10_24_35_btcc)


# Plotting

NEB_plots <- lapply(
  seq(1,3),
  function(i) plot(NEB_eul_list[[i]],
                   quantities = list(fontsize=14), fills = F, labels = F,
                   edges = list(col="royalblue1", lty="solid", lwd=2))
)

SMART10_plots <- lapply(
  seq(1,3),
  function(i) plot(SMART10_eul_list[[i]],
                   quantities = list(fontsize=14), fills = F, labels = F,
                   edges = list(col="maroon2", lty="solid", lwd=2))
)

SMART100_plots <- lapply(
  seq(1,3),
  function(i) plot(SMART100_eul_list[[i]],
                   quantities = list(fontsize=14), fills = F, labels = F,
                   edges = list(col="olivedrab3", lty="solid", lwd=2))
)

PLT_fig3_C <- plot_grid(plotlist = c(NEB_plots, SMART10_plots, SMART100_plots), ncol=3)


# Plotting

figure_3B <- ggplot(DF_forplot$BTBFCC[DF_forplot$BTBFCC$meanCOV.y<9999,], aes(meanCOV.y, AI.y)) +
  geom_point(aes(color=BT.y), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  facet_grid(libprep ~ BT.x) +
  theme_bw() +
  labs(x = "Total Allele Counts", y = "Gene AI - Replicate 2", color = "Binomial Test") +
  scale_x_continuous(trans='log2') +
  theme(legend.position="None", text = element_text(size=18))

figure_3C <- PLT_fig3_C # see inSanityCheck.R to get this figure

figure_3D <- ggplot(CC_DF, aes(y=variable, x=value, col=variable)) +
  geom_jitter() +
  scale_color_manual(values=c("royalblue1","maroon2","olivedrab3")) +
  theme_bw() +
  theme(legend.position = "None") +
  scale_y_discrete(limits = rev(levels(CC_DF$variable))) +
  ylab("Library preparation method") +
  xlab("QCC")

figure_3E <- ggplot(df, aes(x=libprep, y=P, col=libprep)) +
  geom_boxplot() +
  facet_grid(what ~ .) +
  xlab("") +
  ylab("Concordance % \nbetween two replicates") +
  theme_bw() +
  theme(legend.position="bottom", text = element_text(size=18)) +
  guides(col=guide_legend(title="", nrow = 1)) +
  scale_color_manual(labels=methods,
                     values=c("royalblue1","maroon2","olivedrab3")) +
  #scale_color_manual(values=c("#999999", "#66FFB2"), labels=c("no correction", "QCC correction")) +
  scale_x_discrete(labels=c("", "", ""), limits = rev(levels(df$libprep))) +
  coord_flip()


# figure_3D <- ggplot(res62_df_all, aes(x=method, y=FP_rate, col=experiment)) +
#   geom_boxplot() +
#   xlab("") +
#   ylab("False positive rate") +
#   theme_bw() +
#   theme(legend.position="None", text = element_text(size=18)) +
#   scale_color_manual(labels=methods,
#                      values=wes_palette(n=3, name="FantasticFox1")) +
#   scale_x_discrete(labels=c("no correction", "QCC correction"))

PLT_fig3 = plot_grid(
  plot_grid(NULL, figure_3B, figure_3C, labels = c("A", "B", "C"), rel_widths = c(0.3, 0.8, 0.8), nrow = 1),
  plot_grid(figure_3D, figure_3E, labels = c("D", "E"), rel_widths = c(0.8, 0.8)),
  nrow = 2,
  scale = c(0.9, 0.9, 0.7, 0.9, 0.9)
)

cowplot::save_plot(PLT_fig3, file="~/Dropbox (Partners HealthCare)/replicates_ASE/manuscript/Figures/fig.3_v2.pdf", base_height = 10, base_width = 16)
