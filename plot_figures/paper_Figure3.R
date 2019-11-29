setwd("~/Dropbox (Partners HealthCare)/code/ASE/plot_figures")

library(tidyverse)
library(cowplot)
library(eulerr)
library(gridExtra)
library(ggrepel)
library(wesanderson)
library(eulerr)

source("../R/ASE_functions.R")
source("../R/PerformAIAnalysis_CC.v0.3.R")
source("../R/DownstreamASE.R")
source("../plot_figures/utilPlots.R")
source("../plot_figures/takedata30mln.R")

methods = c("NEBNext (100ng)", "SMARTseq (10ng)", "SMARTseq (0.1ng)")

Ntrials = 1
set.seed(1)
#sampleMreps = sapply(1:(max(Ntrials, 10)), function(i){(sample(1:6, 2)-1)*5 + sample(5, 2, replace=T)})
list_of_datas = data_XXmln_noX
list_of_consts = list(CC[[1]], CC[[2]], CC[[3]])
list_of_libprepnames = list("NEBNext (100ng)", "SMARTseq (10ng)", "SMARTseq (0.1ng)")
#reppair = sampleMreps[, 1]
reppair = c(16, 21)
reppair_consts = lapply(1:3, function(i){list_of_consts[[i]][13]})

CC_df <- data.frame(t(rbind(CC[[1]], CC[[2]], CC[[3]])))
colnames(CC_df) <- c("NEBNext (100ng)", "SMARTseq (10ng)", "SMARTseq (0.1ng)")
CC_DF <- reshape2::melt(CC_df)

# -------------------------------------FIG 3 b--------------------------------------------------------

DF_forplot = CreateForplotDF_btNbtcc(lapply(1:3, function(i){
  CreateForplotDF(list_of_datas[[i]], reppair, reppair_consts[[i]], list_of_libprepnames[[i]])
}))

DF_forplot$BTBFCC$BT.x = as.factor(DF_forplot$BTBFCC$BT.x)
levels(DF_forplot$BTBFCC$BT.x) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")
DF_forplot$BTBFCC$BT.y = as.factor(DF_forplot$BTBFCC$BT.y)
levels(DF_forplot$BTBFCC$BT.y) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")

# -------------------------------------FIG 3 d--------------------------------------------------------

# Calculate disagreement percentages

# For before correction

pdis = lapply(1:3, function(m){
  sapply(1:ncol(combn(6, 2)), function(j){
    sample2reps30mln =(combn(6, 2)[, j]-1)*5 + sample(1:5, 2, replace=T)

    sample2reps30mln_6 = sort((sample2reps30mln-1) %/% 5 + 1)
    CCn_sample2reps30mln = sum((5:1)[0:(sample2reps30mln_6[1]-1)]) + sample2reps30mln_6[2]-sample2reps30mln_6[1]

    CC_sample2_30mln = CC[[m]][CCn_sample2reps30mln]

    data_XXmln_noX_sample2 = data_XXmln_noX[[m]][, sort(c(1, sample2reps30mln*2, sample2reps30mln*2+1))]

    p = GetPercDiffForExperiment_BT(df=data_XXmln_noX_sample2, CC_df=CC_sample2_30mln, exp_name = "any")

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

    data_XXmln_noX_sample2 = data_XXmln_noX[[m]][, sort(c(1, sample2reps30mln*2, sample2reps30mln*2+1))]

    p = GetPercDiffForExperiment_BT_CC(df=data_XXmln_noX_sample2, CC_df=CC_sample2_30mln, exp_name="any")

    return(p)
  })
})

df = rbind(
  do.call(rbind, lapply(1:3, function(i){data.frame(Pc=as.vector(pdis[[i]][1,]), Pd=as.vector(pdis[[i]][2:3,]),
                                                    libprep=list_of_libprepnames[[i]], what="binomial")})),
  do.call(rbind, lapply(1:3, function(i){data.frame(Pc=as.vector(pdis_CC[[i]][1,]), Pd=as.vector(pdis_CC[[i]][2:3,]),
                                                    libprep=list_of_libprepnames[[i]], what="corrected")}))
)

df$libprep <- as.factor(df$libprep)
df$libprep <- factor(df$libprep, levels=c("NEBNext (100ng)","SMARTseq (10ng)","SMARTseq (0.1ng)"))

df_concordance <- df


# Figure 3C

res62_df <-
  lapply(1:3, function(m){
    df <- sapply(1:ncol(combn(6, 2)), function(j){
      pairs_sample = (combn(6, 2)[, j]-1)*5 + sample(1:5, 2, replace=T)
      for_6_df <- CountsToAI(data_XXmln_noX[[m]], meth="mergedToProportion")$AI
      res62 <- PerformBinTestAIAnalysisForConditionNPointVect_knownCC(data_XXmln_noX[[m]], vectReps=pairs_sample,
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

# Euler Fig A

CalculateTestWithHalf <- function(df, cc, exp_name, thr=10){
  res = Reduce(function(x, y) merge(x, y, by="ID"),
               lapply(1:(ncol(df)%/%(2*5)), function(j){
                 sample1rep = (j-1)*5 + sample(1:5, 1)
                 PerformBinTestAIAnalysisForConditionNPoint_knownCC(df, sample1rep, cc, thr=thr)[, c("ID", "BT", "BT_CC")]
               })
  )
  names(res) = c("ID",
                 paste0(c("BT.", "BT_CC."),
                        rep(1:(ncol(df)%/%(2*5)), each=2)
                 )
  )
  return(res)
}


# -------------------------------------FIG 3 a--------------------------------------------------------

list_of_constsM = list(mean(CC[[1]][6:15]), mean(CC[[2]]), mean(CC[[3]]))
list_of_libprepnames = list("NEBNext (100ng)", "SMARTseq (10ng)", "SMARTseq (0.1ng)")

df_NEB = CalculateTestWithHalf(data_XXmln_noX[[1]], list_of_constsM[[1]], list_of_libprepnames[[1]], thr=10)
df_SM10 = CalculateTestWithHalf(data_XXmln_noX[[2]], list_of_constsM[[2]], list_of_libprepnames[[2]], thr=10)
df_SM100 = CalculateTestWithHalf(data_XXmln_noX[[3]], list_of_constsM[[3]], list_of_libprepnames[[3]], thr=10)

df_NEB_SM10_SM100 = Reduce(function(x, y) merge(x, y, by="ID"), list(df_NEB, df_SM10, df_SM100))
dft = df_NEB_SM10_SM100
names(dft) = c("ID", paste0("rep", rep(1:18, each=2), '_', c("bt", "btCC")))
rownames(dft) = dft$ID


dft_2 = rbind(dft, dft, dft)
dft_2[1:nrow(dft), c(12+(2:13), 24+(2:13))] = FALSE
dft_2[(nrow(dft)+1):(2*nrow(dft)), c(2:13, 24+(2:13))] = FALSE
dft_2[(nrow(dft)*2+1):(3*nrow(dft)), c(2:13, 12+(2:13))] = FALSE
rownames(dft_2) = c(paste0(rownames(dft), ".1"), paste0(rownames(dft), ".2"), paste0(rownames(dft), ".3"))

eulij3exp = euler(na.omit(dft_2[, c('rep2_bt', 'rep2_btCC', 'rep4_bt', 'rep4_btCC',
                                    'rep8_bt', 'rep8_btCC', 'rep9_bt', 'rep9_btCC',
                                    'rep14_bt', 'rep14_btCC', 'rep15_bt', 'rep15_btCC'
)]),
shape='ellipse')
# eulij3exp_plt = plot(eulij3exp,
#      quantities = F,
#      fills = F,
#      edges = list(col=c("royalblue1","royalblue1","royalblue1","royalblue1",
#                         "maroon2","maroon2","maroon2","maroon2",
#                         "olivedrab3","olivedrab3","olivedrab3","olivedrab3"),
#                   lty=c("dashed","solid","dashed","solid",
#                         "dashed","solid","dashed","solid",
#                         "dashed","solid","dashed","solid"),
#                   lwd=3),
#      labels = F,
#      legend = F)

eulij3exp_copy = eulij3exp

eulij3exp_copy$ellipses[, 'h'] = eulij3exp_copy$ellipses[, 'h'] - min(eulij3exp_copy$ellipses[, 'h'] - eulij3exp_copy$ellipses[, 'a']) + 1
eulij3exp_copy$ellipses[, 'k'] = eulij3exp_copy$ellipses[, 'k'] - min(eulij3exp_copy$ellipses[, 'k'] - eulij3exp_copy$ellipses[, 'a']) + 1

circles_b = data.frame(
  left_b  = eulij3exp_copy$ellipses[, 'h'] - eulij3exp_copy$ellipses[, 'a'],
  right_b = eulij3exp_copy$ellipses[, 'h'] + eulij3exp_copy$ellipses[, 'a'],
  lower_b = eulij3exp_copy$ellipses[, 'k'] - eulij3exp_copy$ellipses[, 'a'],
  upper_b = eulij3exp_copy$ellipses[, 'k'] + eulij3exp_copy$ellipses[, 'a']
)
euler_b = data.frame(
  left_b  = c(min(circles_b$left_b[1:4]), min(circles_b$left_b[5:8]), min(circles_b$left_b[9:12])),
  right_b = c(max(circles_b$right_b[1:4]), max(circles_b$right_b[5:8]), max(circles_b$right_b[9:12])),
  lower_b = c(min(circles_b$lower_b[1:4]), min(circles_b$lower_b[5:8]), min(circles_b$lower_b[9:12])),
  upper_b = c(max(circles_b$upper_b[1:4]), max(circles_b$upper_b[5:8]), max(circles_b$upper_b[9:12]))
)
euler_c = t(sapply(1:3, function(i){c(mean(as.numeric(euler_b[i, 1:2])), mean(as.numeric(euler_b[i, 3:4])))}))
euler_s = sapply(1:3, function(i){max(abs(euler_b[i, 1]-euler_b[i, 2]), abs(euler_b[i, 3]-euler_b[i, 4]))})


Meuler_s = max(euler_s)
euler_c_new = cbind(c(Meuler_s/2, Meuler_s/2, Meuler_s/2),
                    c(5*Meuler_s/2, 3*Meuler_s/2, Meuler_s/2))
euler_c_shift = euler_c_new - euler_c

euler_c_shift_x4 = do.call(rbind, lapply(1:3, function(i){
  rbind(euler_c_shift[i,],euler_c_shift[i,],euler_c_shift[i,],euler_c_shift[i,])
}))

eulij3exp_copy$ellipses[, c('h','k')] = eulij3exp_copy$ellipses[, c('h','k')] + euler_c_shift_x4

eulij3exp_plt = plot(eulij3exp_copy,
     quantities = F,
     fills = F,
     edges = list(col=c("royalblue1","royalblue1","royalblue1","royalblue1",
                        "maroon2","maroon2","maroon2","maroon2",
                        "olivedrab3","olivedrab3","olivedrab3","olivedrab3"),
                  lty=c("dashed","solid","dashed","solid",
                        "dashed","solid","dashed","solid",
                        "dashed","solid","dashed","solid"),
                  lwd=3,
                  alpha = c(0.5,1,0.5,1)),
     labels = F,
     legend = F)
# -------------------------------------FIG 3 e--------------------------------------------------------

# Data prep

#NEB

## single rep
df_NEB_2_3 = CalculateTestWithHalf1(data_XXmln_noX[[1]][, c(1, 2*5*2, 2*5*2+1, 3*5*2, 3*5*2+1)],
                                    CC[[1]][6], list_of_libprepnames[[1]], thr=10)
## merge two reps
df_NEB_24 = CalculateTestWithHalf1(MergeSumCounts(data_XXmln_noX[[1]][, c(1, 2*5*2, 2*5*2+1, 4*5*2, 4*5*2+1)]),
                                   CC[[1]][7], list_of_libprepnames[[1]], thr=10)
df_NEB_35 = CalculateTestWithHalf1(MergeSumCounts(data_XXmln_noX[[1]][, c(1, 3*5*2, 3*5*2+1, 5*5*2, 5*5*2+1)]),
                                   CC[[1]][11], list_of_libprepnames[[1]], thr=10)
df_NEB_24_35 = merge(df_NEB_24, df_NEB_35, by="ID")


#SMART100

## single rep
df_SM100_2_3 = CalculateTestWithHalf1(data_XXmln_noX[[3]][, c(1, 2*5*2, 2*5*2+1, 3*5*2, 3*5*2+1)],
                                      CC[[3]][6], list_of_libprepnames[[3]], thr=10)
## merge two reps
df_SM100_24 = CalculateTestWithHalf1(MergeSumCounts(data_XXmln_noX[[3]][, c(1, 2*5*2, 2*5*2+1, 4*5*2, 4*5*2+1)]),
                                     CC[[3]][7], list_of_libprepnames[[3]], thr=10)
df_SM100_35 = CalculateTestWithHalf1(MergeSumCounts(data_XXmln_noX[[3]][, c(1, 3*5*2, 3*5*2+1, 5*5*2, 5*5*2+1)]),
                                     CC[[3]][11], list_of_libprepnames[[3]], thr=10)
df_SM100_24_35 = merge(df_SM100_24, df_SM100_35, by="ID")


#SMART10

## single rep
df_SM10_2_3 = CalculateTestWithHalf1(data_XXmln_noX[[2]][, c(1, 2*5*2, 2*5*2+1, 3*5*2, 3*5*2+1)],
                                     CC[[2]][6], list_of_libprepnames[[2]], thr=10)
## merge two reps
df_SM10_24 = CalculateTestWithHalf1(MergeSumCounts(data_XXmln_noX[[2]][, c(1, 2*5*2, 2*5*2+1, 4*5*2, 4*5*2+1)]),
                                    CC[[2]][7], list_of_libprepnames[[2]], thr=10)
df_SM10_35 = CalculateTestWithHalf1(MergeSumCounts(data_XXmln_noX[[2]][, c(1, 3*5*2, 3*5*2+1, 5*5*2, 5*5*2+1)]),
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


# # -------------------------------------FIG 3 b--------------------------------------------------------
# data_30mln_noX_sample2 = lapply(data_XXmln_noX, function(x){
#   x[, sort(c(1, reppair*2, reppair*2+1))]
# })
#
# DF_BTCC = do.call(rbind, lapply(1:3, function(i){
#   GetDataForExperiment_BTCC(df = data_30mln_noX_sample2, CC_df = reppair_consts,
#                             exp_name = list_of_libprepnames[[i]], exp_n = i)
# }))
#
# ggplot(DF_BTCC[DF_BTCC$meanCOV.y < 9999, ], aes(meanCOV.y, AI.y)) +
#   geom_point(aes(color=BT.y), size=0.8) +
#   scale_color_manual(values=c("salmon","royalblue1")) +
#   facet_grid(libprep ~ BT.x) +
#   theme_bw() +
#   labs(x = "Total Allelic Counts", y = "Gene AI - Replicate 2", color = "Binomial Test") +
#   scale_x_continuous(trans='log10') +
#   theme(legend.position="None", text = element_text(size=Tsize),
#         strip.background = element_rect(fill="#E0E0E0"))

#===========================================================================================================
# Plotting

NEB_plots <- lapply(
  seq(1,3),
  function(i) plot(NEB_eul_list[[i]],
                   quantities = list(fontsize=14), fills = F, labels = F,
                   edges = list(col="royalblue1", lty="solid", lwd=3))
)

SMART10_plots <- lapply(
  seq(1,3),
  function(i) plot(SMART10_eul_list[[i]],
                   quantities = list(fontsize=14), fills = F, labels = F,
                   edges = list(col="maroon2", lty="solid", lwd=3))
)

SMART100_plots <- lapply(
  seq(1,3),
  function(i) plot(SMART100_eul_list[[i]],
                   quantities = list(fontsize=14), fills = F, labels = F,
                   edges = list(col="olivedrab3", lty="solid", lwd=3))
)

figure_3E <- plot_grid(plotlist = c(NEB_plots, SMART10_plots, SMART100_plots), ncol=3)
                        #labels = c("A", "B", "C"),
                        #label_x = .5, hjust = 0)

Tsize = 18

figure_3A <- eulij3exp_plt

figure_3B <- ggplot(DF_forplot$BTBFCC[DF_forplot$BTBFCC$meanCOV.y<9999,], aes(meanCOV.y, AI.y)) +
  geom_point(aes(color=BT.y), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  scale_y_continuous(breaks = pretty(DF_forplot$BTBFCC[DF_forplot$BTBFCC$meanCOV.y<9999,]$AI.y, n = 3)) +
  facet_grid(libprep ~ BT.x) +
  theme_bw() +
  labs(x = "Total Allele Counts", y = "Gene AI - Replicate 2", color = "Binomial Test") +
  scale_x_continuous(trans='log2') +
  theme(legend.position="None", text = element_text(size=Tsize), strip.text.x = element_text(size=14), strip.text.y = element_text(size=14))

figure_3C <- ggplot(CC_DF, aes(y=variable, x=value, col=variable)) +
  geom_jitter() +
  scale_color_manual(values=c("royalblue1","maroon2","olivedrab3")) +
  theme_bw() +
  xlim(c(1,3)) +
  theme(legend.position = "None", text = element_text(size=Tsize)) +
  scale_y_discrete(limits = rev(levels(CC_DF$variable)), labels = c("", "", "")) +
  ylab("Experiment") +
  xlab("QCC")

figure_3D <- ggplot(df_concordance, aes(x=libprep, y=Pc, col=libprep)) +
  geom_boxplot() +
  facet_grid(what ~ .) +
  xlab("Experiment") +
  ylab("Concordance between two replicates %") +
  theme_bw() +
  theme(legend.position="None", text = element_text(size=Tsize)) +
  scale_color_manual(labels=methods,
                     values=c("royalblue1","maroon2","olivedrab3")) +
  scale_x_discrete(labels=c("", "", ""), limits = rev(levels(df_concordance$libprep))) +
  coord_flip()

PLT_fig3 = plot_grid(
  plot_grid(
    plot_grid(figure_3A, figure_3C, figure_3D, labels = c("A", "C", "D"), rel_heights = c(0.45, 0.22, 0.33), scale = c(0.9, 0.9, 0.9), ncol = 1),
    plot_grid(figure_3B, labels = c("B"), scale = c(0.97)),
    ncol=2
  ),
  plot_grid(figure_3E, labels = c("E"), scale = c(0.8)),
  rel_heights=c(63,37), ncol=1
)

cowplot::save_plot(PLT_fig3, file="fig.3_v5.pdf", base_height = 20, base_width = 12)

# numbers for fig3a :
df = cbind(
  NEB_100ng_bt = euler(na.omit(dft_2[, c('rep2_bt', 'rep4_bt')]), shape='ellipse')$original,
  NEB_100ng_bt_QCC = euler(na.omit(dft_2[, c('rep2_btCC', 'rep4_btCC')]), shape='ellipse')$original,
  SMART_10ng_bt = euler(na.omit(dft_2[, c('rep8_bt', 'rep9_bt')]), shape='ellipse')$original,
  SMART_10ng_bt_QCC = euler(na.omit(dft_2[, c('rep8_btCC', 'rep9_btCC')]), shape='ellipse')$original,
  SMART_0.1ng_bt = euler(na.omit(dft_2[, c('rep14_bt', 'rep15_bt')]), shape='ellipse')$original,
  SMART_0.1ng_bt_QCC = euler(na.omit(dft_2[, c('rep14_btCC', 'rep15_btCC')]), shape='ellipse')$original
)
row.names(df) = c("only1", "only2", "intersection")
df_NP = data.frame(N = df[3,], P = paste(round(df[3,]/(df[1,]+df[2,]+df[3,]) * 100, 1), "%"))
plot(tableGrob(df_NP, theme = ttheme_minimal()))


PLT_fig3 = plot_grid(
  plot_grid(
    plot_grid(plot_grid(figure_3A, tableGrob(df_NP, theme = ttheme_minimal()), ncol=2),
              figure_3C, figure_3D, labels = c("A", "C", "D"), rel_heights = c(0.45, 0.22, 0.33), scale = c(0.9, 0.9, 0.9), ncol = 1),
    plot_grid(figure_3B, labels = c("B"), scale = c(0.97)),
    ncol=2
  ),
  plot_grid(figure_3E, labels = c("E"), scale = c(0.8)),
  rel_heights=c(63,37), ncol=1
)

cowplot::save_plot(PLT_fig3, file="fig.3_v5.pdf", base_height = 20, base_width = 12)

# PLT_fig3 = plot_grid(
#   plot_grid(figure_3A, figure_3C, labels = c("A", "B"), rel_widths = c(0.4, 0.8), nrow = 1, scale = c(0.9, 0.7)),
#   plot_grid(figure_3D, figure_3E, labels = c("D", "E"), rel_widths = c(0.8, 0.8)),
#   nrow = 2,
#   scale = c(0.9, 0.9, 0.9, 0.9, 0.9)
# )
#
# cowplot::save_plot(PLT_fig3, file="~/Dropbox (Partners HealthCare)/replicates_ASE/manuscript/Figures/fig.3_v4.pdf", base_height = 10, base_width = 16)
# cowplot::save_plot(figure_3B, file="~/Dropbox (Partners HealthCare)/replicates_ASE/manuscript/Figures/fig.3C_.pdf", base_height = 8, base_width = 4)




# concordance-discordance:

# figure_3E_c <- ggplot(df_concordance, aes(x=libprep, y=Pc, col=libprep)) +
#   geom_boxplot() + geom_point() +
#   facet_grid(~what) +
#   xlab("") +
#   ylab("Concordance % \nbetween two replicates") +
#   theme_bw() +
#   theme(legend.position="right", text = element_text(size=Tsize)) +
#   guides(col=guide_legend(title="")) +
#   scale_color_manual(labels=methods,
#                      values=c("royalblue1","maroon2","olivedrab3")) +
#   #scale_color_manual(values=c("#999999", "#66FFB2"), labels=c("no correction", "QCC correction")) +
#   scale_x_discrete(labels=c("", "", ""), limits = (levels(df_concordance$libprep))) #+
#   #coord_flip()
#
# figure_3E_d <- ggplot(df_concordance, aes(x=libprep, y=Pd, col=libprep)) +
#   geom_boxplot() + geom_point() +
#   facet_grid(~what) +
#   xlab("") +
#   ylab("Discovery Desagreement Rate \nbetween two replicates") +
#   theme_bw() +
#   theme(legend.position="bottom", text = element_text(size=Tsize)) +
#   guides(col=guide_legend(title="", nrow = 1)) +
#   scale_color_manual(labels=methods,
#                      values=c("royalblue1","maroon2","olivedrab3")) +
#   #scale_color_manual(values=c("#999999", "#66FFB2"), labels=c("no correction", "QCC correction")) +
#   scale_x_discrete(labels=c("", "", ""), limits = (levels(df_concordance$libprep))) #+
#   #coord_flip()
#
# cowplot::plot_grid(figure_3E_d, figure_3E_c)

# figure_3D <- ggplot(res62_df_all, aes(x=method, y=FP_rate, col=experiment)) +
#   geom_boxplot() +
#   xlab("") +
#   ylab("False positive rate") +
#   theme_bw() +
#   theme(legend.position="None", text = element_text(size=Tsize)) +
#   scale_color_manual(labels=methods,
#                      values=wes_palette(n=3, name="FantasticFox1")) +
#   scale_x_discrete(labels=c("no correction", "QCC correction"))
