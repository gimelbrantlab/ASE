setwd("~/Dropbox (Partners HealthCare)/code/ASE/plot_figures")

library(tidyverse)
library(fitdistrplus)
library(cowplot)
library(ggrepel)
library(wesanderson)
library(gridExtra)
library(gplots)

#library(DESeq2)
library(Biobase)

source("../R/ASE_functions.R")
source("../R/PerformAIAnalysis_CC.R")
source("../R/DownstreamASE.R")
source("../plot_figures/utilPlots.R")

source("../../../variation_project/GTEx/GTExR/R/getAI_Data.R")
source("../../../variation_project/GTEx/GTExR/R/utils.R")


tissues  = c("Kidney - Cortex", "Heart - Left Ventricle", "Liver", "Lung", "Pancreas", "Stomach", "Small Intestine - Terminal Ileum", "Spleen")

all_tissues = lapply(tissues, function(x){get_ASE_files_list(x)})
names(all_tissues) = tissues

merged_tissues = Reduce(function(x,y){merge(x,y,by="id",all=T)},
       lapply(1:length(tissues), function(i){
         df = data.frame(id=all_tissues[[i]]$id, tissue = 1)
         names(df) = c("id", tissues[i])
         return(df)
       }
       ))


## "GTEX-11NUK"

s1_ll = "GTEX-11NUK"
tab_liver = read_counts("Liver")
tab_lung = read_counts("Lung")

tab_liver_s1 = tab_liver[startsWith(tab_liver$sample_id, s1_ll), ]
tab_lung_s1 = tab_lung[startsWith(tab_lung$sample_id, s1_ll), ]

print(unique(tab_liver_s1$sample_id))
print(unique(tab_lung_s1$sample_id))


tab_s1_ll = merge(tab_liver_s1[,1:3], tab_lung_s1[,1:3], by="GENE_ID")
names(tab_s1_ll)[1] = "ensembl_gene_id"


simCC_GTEx = seq(1,3,0.05)

interrep_CCsim_GTEX_diff = lapply(simCC_GTEx, function(cc){
  PerformBinTestAIAnalysisForTwoConditions_knownCC(tab_s1_ll,
                                                   vect1CondReps=1,
                                                   vect2CondReps=2,
                                                   vect1CondRepsCombsCC=cc, vect2CondRepsCombsCC=cc,
                                                   Q=0.95, thr=8, thrUP=NA, thrType="each", minDifference=NA)
})

interrep_CCsim_GTEX_diff_df = do.call(rbind, lapply(1:length(interrep_CCsim_GTEX_diff), function(x){
  df = na.omit(interrep_CCsim_GTEX_diff[[x]][, c("BT","BT_CC")])
  data.frame(
    DiffASE_proportion = sum(df[, 2])/nrow(df),
    DiffASE_N = sum(df[, 2]),
    gene_num = nrow(df),
    QCC = simCC_GTEx[x]
  )
}))

ggF4E = ggplot(interrep_CCsim_GTEX_diff_df,
               aes(QCC, DiffASE_N, group=QCC)) +
  geom_point() +
  xlab("Tested QCC") +
  ylab(paste("Number of genes \ncalled differential /", interrep_CCsim_GTEX_diff_df$gene_num[1])) +
  theme_bw() +
  scale_y_continuous(breaks=seq(0, max(interrep_CCsim_GTEX_diff_df$DiffASE_N)+1, 2),
                     limits = c(0, max(interrep_CCsim_GTEX_diff_df$DiffASE_N)+1)) +
  theme(legend.position="None", text = element_text(size=18))




df_testhalf_1 = lapply(simCC_GTEx, function(cc){
  PerformBinTestAIAnalysisForConditionNPoint_knownCC(tab_s1_ll, vectReps = 1,
                                                     vectRepsCombsCC = cc,
                                                     pt = 0.5, binNObs=40, Q=0.95,
                                                     thr=8, thrUP=NA, thrType="each", minDifference=NA)
})

df_unbiased_1 = rbind(
  data.frame(df_testhalf_1[[1]][df_testhalf_1[[1]]$BT_CC == F, ], class = "Unbiased", tissue = "Liver")
)
df_biased_1 = do.call(rbind, lapply(seq(1,31,4), function(i){
  data.frame(df_testhalf_1[[i]][df_testhalf_1[[i]]$BT_CC == T, ], class = simCC_GTEx[i], tissue = "Liver")
}))

df_testhalf_2 = lapply(simCC_GTEx, function(cc){
  PerformBinTestAIAnalysisForConditionNPoint_knownCC(tab_s1_ll, vectReps = 2,
                                                     vectRepsCombsCC = cc,
                                                     pt = 0.5, binNObs=40, Q=0.95,
                                                     thr=8, thrUP=NA, thrType="each", minDifference=NA)
})

df_unbiased_2 = rbind(
  data.frame(df_testhalf_2[[1]][df_testhalf_2[[1]]$BT_CC == F, ], class = "Unbiased", tissue = "Lung")
)
df_biased_2 = do.call(rbind, lapply(seq(1,31,4), function(i){
  data.frame(df_testhalf_2[[i]][df_testhalf_2[[i]]$BT_CC == T, ], class = simCC_GTEx[i], tissue = "Lung")
}))

ggF3C = ggplot() +
  geom_point(data=rbind(df_unbiased_1, df_unbiased_2), aes(sumCOV, AI), col='grey', size=0.8) +
  geom_point(data=rbind(df_biased_1, df_biased_2),   aes(sumCOV, AI, col=class),  size=0.8) +
  facet_grid(~ tissue) +
  theme_bw() +
  scale_x_log10() +
  xlab("AI") + ylab("Total Gene Coverage") + labs(col = "Biased:\nQCC ") +
  theme(legend.position="right", text=element_text(size=18))

SuppFig_GTEx_simCC = plot_grid(ggF3C, ggF4E, ncol=1, labels=c("C", "D"), rel_heights = c(1,0.7))
SuppFig_GTEx_simCC

# ggplot(df_unbiased_1, aes(sumCOV)) +
#   geom_histogram(binwidth = 10) +
#   theme_bw() + xlim(0, 1000) +
#   geom_vline(xintercept = quantile(df_unbiased_1$sumCOV, c(0.5,0.75,0.95)))


cowplot::save_plot(SuppFig_GTEx_simCC, file="SuppFig_GTEx.pdf", base_height = 8, base_width = 6)








# LUNG and LIVER:

samples_lung_liver = na.omit(merged_tissues[, c("id", 'Lung', "Liver")])


ss_ll = c(1:5, sample(6:nrow(samples_lung_liver), 5))
samples_lung_liver$id[ss_ll]

lung_liver_smpl5_list = lapply(1:5, function(i){

  s1_ll = samples_lung_liver$id[i]


  tab_liver = read_counts("Liver")
  tab_lung = read_counts("Lung")

  tab_liver_s1 = tab_liver[startsWith(tab_liver$sample_id, s1_ll), ]
  tab_lung_s1 = tab_lung[startsWith(tab_lung$sample_id, s1_ll), ]

  print(unique(tab_liver_s1$sample_id))
  print(unique(tab_lung_s1$sample_id))


  tab_s1_ll = merge(tab_liver_s1[,1:3], tab_lung_s1[,1:3], by="GENE_ID")
  names(tab_s1_ll)[1] = "ensembl_gene_id"


  simCC_GTEx = seq(1,3,0.05)

  interrep_CCsim_GTEX_diff = lapply(simCC_GTEx, function(cc){
    PerformBinTestAIAnalysisForTwoConditions_knownCC(tab_s1_ll,
                                                     vect1CondReps=1,
                                                     vect2CondReps=2,
                                                     vect1CondRepsCombsCC=cc, vect2CondRepsCombsCC=cc,
                                                     Q=0.95, thr=8, thrUP=NA, thrType="each", minDifference=NA)
  })

  interrep_CCsim_GTEX_diff_df = do.call(rbind, lapply(1:length(interrep_CCsim_GTEX_diff), function(x){
    df = na.omit(interrep_CCsim_GTEX_diff[[x]][, c("BT","BT_CC")])
    data.frame(
      DiffASE_proportion = sum(df[, 2])/nrow(df),
      DiffASE_N = sum(df[, 2]),
      gene_num = nrow(df),
      QCC = simCC_GTEx[x]
    )
  }))

  ggF4E = ggplot(interrep_CCsim_GTEX_diff_df,
         aes(QCC, DiffASE_N, group=QCC)) +
    geom_point() +
    xlab("Tested QCC") +
    ylab(paste("Number of genes \ncalled differential / ", interrep_CCsim_GTEX_diff_df$gene_num[1])) +
    theme_bw() +
    scale_y_continuous(breaks=seq(0, max(interrep_CCsim_GTEX_diff_df$DiffASE_N)+1, 1)) +
    theme(legend.position="None", text = element_text(size=18))


  interrep_CCsim_GTEX_diff_all = do.call(rbind, lapply(1:length(interrep_CCsim_GTEX_diff),
                                                       function(i){data.frame(interrep_CCsim_GTEX_diff[[i]], QCC = simCC_GTEx[i])}))



  ggDiff = plot_grid(
    ggplot(interrep_CCsim_GTEX_diff_all, aes(sumCOV_1, AI_1)) +
      geom_point(col='grey', size=0.5) +
      geom_point(data=interrep_CCsim_GTEX_diff_all[interrep_CCsim_GTEX_diff_all$BT_CC == T, ],
                 aes(sumCOV_1, AI_1), col='red', size=1) +
      xlab("Coverage Liver") +
      ylab("AI Liver") +
      theme_bw() +
      scale_x_log10() +
      facet_grid(~ QCC) +
      theme(legend.position="None", text = element_text(size=10))
    ,
    ggplot(interrep_CCsim_GTEX_diff_all, aes(sumCOV_1, AI_2)) +
      geom_point(col='grey', size=0.5) +
      geom_point(data=interrep_CCsim_GTEX_diff_all[interrep_CCsim_GTEX_diff_all$BT_CC == T, ],
                 aes(sumCOV_1, AI_2), col='red', size=1) +
      xlab("Coverage Lung") +
      ylab("AI Liver") +
      theme_bw() +
      scale_x_log10() +
      facet_grid(~ QCC) +
      theme(legend.position="None", text = element_text(size=10))
    ,
    ggplot(interrep_CCsim_GTEX_diff_all, aes(sumCOV_1, AI_1-AI_2)) +
      geom_point(col='grey', size=0.5) +
      geom_point(data=interrep_CCsim_GTEX_diff_all[interrep_CCsim_GTEX_diff_all$BT_CC == T, ],
                 aes(sumCOV_1, AI_1-AI_2), col='red', size=1) +
      xlab("Coverage Liver") +
      ylab("AI Liver - AI Lung") +
      theme_bw() +
      scale_x_log10() +
      facet_grid(~ QCC) +
      theme(legend.position="None", text = element_text(size=10))
    ,
    ncol=1
  )

  return(list(ggF4E, ggDiff, interrep_CCsim_GTEX_diff_df, interrep_CCsim_GTEX_diff_all))
})
names(lung_liver_smpl5_list) = samples_lung_liver$id[ss_ll][1:5]


interrep_CCsim_GTEX_diff_df_X5 = do.call(rbind, lapply(1:5, function(i){
  data.frame(lung_liver_smpl5_list[[i]][3], ID = names(lung_liver_smpl5_list)[i])
}))
gg_fig4E_ll = ggplot(interrep_CCsim_GTEX_diff_df_X5,
       aes(QCC, DiffASE_N, group=QCC)) +
  geom_point() +
  xlab("Tested QCC") +
  ylab(paste("Number of genes \ncalled differential / ", interrep_CCsim_GTEX_diff_df$gene_num[1])) +
  theme_bw() +
  facet_grid(~ ID) +
  theme(legend.position="None", text = element_text(size=18))




# LUNG and LIVER:

# c("Kidney - Cortex", "Heart - Left Ventricle", "Liver", "Lung", "Pancreas", "Stomach", "Small Intestine - Terminal Ileum", "Spleen")
samples_heart_spleen = na.omit(merged_tissues[, c("id", "Heart - Left Ventricle", "Stomach")])


ss_hs = c(sample(6:nrow(samples_heart_spleen), 5))
samples_heart_spleen$id[ss_hs]

heart_spleen_smpl5_list = lapply(1:5, function(i){

  s1_hs = samples_heart_spleen$id[i]


  tab_heart = read_counts("Heart - Left Ventricle")
  tab_spleen = read_counts("Stomach")

  tab_heart_s1 = tab_heart[startsWith(tab_heart$sample_id, s1_hs), ]
  tab_spleen_s1 = tab_spleen[startsWith(tab_spleen$sample_id, s1_hs), ]

  print(unique(tab_heart_s1$sample_id))
  print(unique(tab_spleen_s1$sample_id))


  tab_s1_hs = merge(tab_heart_s1[,1:3], tab_spleen_s1[,1:3], by="GENE_ID")
  names(tab_s1_hs)[1] = "ensembl_gene_id"


  simCC_GTEx = seq(1,1.6,0.05)

  interrep_CCsim_GTEX_diff = lapply(simCC_GTEx, function(cc){
    PerformBinTestAIAnalysisForTwoConditions_knownCC(tab_s1_hs,
                                                     vect1CondReps=1,
                                                     vect2CondReps=2,
                                                     vect1CondRepsCombsCC=cc, vect2CondRepsCombsCC=cc,
                                                     Q=0.95, thr=8, thrUP=NA, thrType="each", minDifference=NA)
  })

  interrep_CCsim_GTEX_diff_df = do.call(rbind, lapply(1:length(interrep_CCsim_GTEX_diff), function(x){
    df = na.omit(interrep_CCsim_GTEX_diff[[x]][, c("BT","BT_CC")])
    data.frame(
      DiffASE_proportion = sum(df[, 2])/nrow(df),
      DiffASE_N = sum(df[, 2]),
      gene_num = nrow(df),
      QCC = simCC_GTEx[x]
    )
  }))

  ggF4E = ggplot(interrep_CCsim_GTEX_diff_df,
                 aes(QCC, DiffASE_N, group=QCC)) +
    geom_point() +
    xlab("Tested QCC") +
    ylab(paste("Number of genes \ncalled differential / ", interrep_CCsim_GTEX_diff_df$gene_num[1])) +
    theme_bw() +
    theme(legend.position="None", text = element_text(size=18))


  interrep_CCsim_GTEX_diff_all = do.call(rbind, lapply(1:length(interrep_CCsim_GTEX_diff),
                                                       function(i){data.frame(interrep_CCsim_GTEX_diff[[i]], QCC = simCC_GTEx[i])}))



  ggDiff = plot_grid(
    ggplot(interrep_CCsim_GTEX_diff_all, aes(sumCOV_1, AI_1)) +
      geom_point(col='grey', size=0.5) +
      geom_point(data=interrep_CCsim_GTEX_diff_all[interrep_CCsim_GTEX_diff_all$BT_CC == T, ],
                 aes(sumCOV_1, AI_1), col='red', size=1) +
      xlab("Coverage Heart") +
      ylab("AI Heart") +
      theme_bw() +
      scale_x_log10() +
      facet_grid(~ QCC) +
      theme(legend.position="None", text = element_text(size=12))
    ,
    ggplot(interrep_CCsim_GTEX_diff_all, aes(sumCOV_1, AI_2)) +
      geom_point(col='grey', size=0.5) +
      geom_point(data=interrep_CCsim_GTEX_diff_all[interrep_CCsim_GTEX_diff_all$BT_CC == T, ],
                 aes(sumCOV_1, AI_2), col='red', size=1) +
      xlab("Coverage Spleen") +
      ylab("AI Heart") +
      theme_bw() +
      scale_x_log10() +
      facet_grid(~ QCC) +
      theme(legend.position="None", text = element_text(size=12))
    ,
    ggplot(interrep_CCsim_GTEX_diff_all, aes(sumCOV_1, AI_1-AI_2)) +
      geom_point(col='grey', size=0.5) +
      geom_point(data=interrep_CCsim_GTEX_diff_all[interrep_CCsim_GTEX_diff_all$BT_CC == T, ],
                 aes(sumCOV_1, AI_1-AI_2), col='red', size=1) +
      xlab("Coverage Heart") +
      ylab("AI Heart - AI Spleen") +
      theme_bw() +
      scale_x_log10() +
      facet_grid(~ QCC) +
      theme(legend.position="None", text = element_text(size=12))
    ,
    ncol=1
  )

  return(list(ggF4E, ggDiff, interrep_CCsim_GTEX_diff_df, interrep_CCsim_GTEX_diff_all))
})
names(heart_spleen_smpl5_list) = samples_heart_spleen$id[ss_hs][1:5]



interrep_CCsim_GTEX_diff_df_X5_2 = do.call(rbind, lapply(1:5, function(i){
  data.frame(heart_spleen_smpl5_list[[i]][3], ID = names(heart_spleen_smpl5_list)[i])
}))
gg_fig4E_hs = ggplot(interrep_CCsim_GTEX_diff_df_X5_2,
                     aes(QCC, DiffASE_N, group=QCC)) +
  geom_point() +
  xlab("Tested QCC") +
  ylab(paste("Number of genes \ncalled differential / ", interrep_CCsim_GTEX_diff_df$gene_num[1])) +
  theme_bw() +
  facet_grid(~ ID) +
  theme(legend.position="None", text = element_text(size=18))




plot_grid(gg_fig4E_ll, gg_fig4E_hs, nrow =2)



ggplot(interrep_CCsim_GTEX_diff[[1]], aes(sumCOV_1)) +
  geom_histogram(binwidth = 10) + xlim(0,1000)
