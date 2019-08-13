setwd("~/Dropbox (Partners HealthCare)/replicates_ASE/code/ASE/plot_figures")

library(tidyverse)
library(cowplot)
library(eulerr)
library(gridExtra)

source("../R/ASE_functions.R")
source("../R/PerformAIAnalysis_CC.R")
source("../R/DownstreamASE.R")
source("../plot_figures/utilPlots.R")
source("../plot_figures/takedata30mln.R")

reps_6_from30 = list(c(3,8,13,19,24,28),
                     c(4,10,12,19,24,28),
                     c(5,10,12,16,23,27))
# reps_6_from30 = list((1:5)*5+1, (0:5)*5+1, (0:5)*5+1)
list_of_datas = data_30mln_noX; list_of_datas[[1]] = data_30mln_noX[[1]][, -c(2:11)]
list_of_consts = list(CC[[1]], CC[[2]], CC[[3]])
list_of_constsM = list(mean(CC[[1]][6:15]), mean(CC[[2]]), mean(CC[[3]]))
list_of_libprepnames = list("NEBNext (100ng)", "SMARTseq (10ng)", "SMARTseq (0.1ng)")
list_of_Nreps = list(6,6,6)

res62_df <-
  lapply(1:3, function(m){
    for_6_df <- CountsToAI(data_30mln_noX[[m]], reps_6_from30[[m]], meth="mergedToProportion")$AI
    df <- sapply(1:length(list_of_consts[[m]]), function(j){
      pairs_sample = reps_6_from30[[m]][combn(list_of_Nreps[[m]], 2)[, j]]
      res62 <- PerformBinTestAIAnalysisForConditionNPointVect_knownCC(data_30mln_noX[[m]],
                                                                      vectReps=pairs_sample,
                                                                      vectRepsCombsCC = list_of_consts[[m]][j],
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
    df$c2_FP_BTCC_rate <- df$FP_BTCC / df$all_BTCC
    df$a2_FP_BT_rate <- df$FP_BT / df$all_BT
    print(df)
    return(df)
  })

res62_df_NEB <- res62_df[[1]]
res62_df_NEB$experiment <- "NEBNext (100ng)"

res62_df_SMART10 <- res62_df[[2]]
res62_df_SMART10$experiment <- "SMARTseq (10ng)"

res62_df_SMART100 <- res62_df[[3]]
res62_df_SMART100$experiment <- "SMARTseq (0.1ng)"

res62_df_all <- data.frame(rbind(res62_df_NEB, res62_df_SMART10, res62_df_SMART100))
res62_df_all <- res62_df_all[,c(7,5:6)] %>% gather(key="method", value = "FP_rate", -experiment)

res61_df <-
  lapply(1:3, function(m){
    for_6_df <- CountsToAI(data_30mln_noX[[m]], reps_6_from30[[m]], meth="mergedToProportion")$AI
    df <- sapply(1:list_of_Nreps[[m]], function(j){
      res61 <- PerformBinTestAIAnalysisForConditionNPointVect_knownCC(data_30mln_noX[[m]],  #list_of_datas[[m]],
                                                                      vectReps=reps_6_from30[[m]][j],
                                                                      vectRepsCombsCC = list_of_consts[[m]][j],
                                                                      ptVect = for_6_df,
                                                                      thr=10)
      FP_BT1 <- sum(res61$BT, na.rm = T)
      all_BT1 <- length(res61$BT) - sum(is.na(res61$BT))
      return(c(FP_BT1, all_BT1))
    })
    df <- data.frame(t(df))
    colnames(df) <- c("FP_BT1", "all_BT1")
    df$a1_FP_BT1_rate <- df$FP_BT1 / df$all_BT1
    print(df)
    return(df)
  })

res61_df_NEB <- res61_df[[1]]
res61_df_NEB$experiment <- "NEBNext (100ng)"

res61_df_SMART10 <- res61_df[[2]]
res61_df_SMART10$experiment <- "SMARTseq (10ng)"

res61_df_SMART100 <- res61_df[[3]]
res61_df_SMART100$experiment <- "SMARTseq (0.1ng)"

res61_df_all <- data.frame(rbind(res61_df_NEB, res61_df_SMART10, res61_df_SMART100))
res61_df_all <- res61_df_all[,c(4,3)] %>% gather(key="method", value = "FP_rate", -experiment)



###

# Data for Figure 4d

data1c = sapply(reps_6_from30[[1]][2:6], function(i){rowSums(data_30mln_noX[[1]][, c(i*2-1,i*2)])})
data2c = sapply(reps_6_from30[[2]][1:6], function(i){rowSums(data_30mln_noX[[2]][, c(i*2-1,i*2)])})
data3c = sapply(reps_6_from30[[3]][1:6], function(i){rowSums(data_30mln_noX[[3]][, c(i*2-1,i*2)])})


dataC = data.frame(id = data_30mln_noX[[1]]$ensembl_gene_id,
                   mc1 = rowMeans(data1c), sdc1 = apply(data1c, 1, sd), vc1 = apply(data1c, 1, var),
                   mc2 = rowMeans(data2c), sdc2 = apply(data2c, 1, sd), vc2 = apply(data2c, 1, var),
                   mc3 = rowMeans(data3c), sdc3 = apply(data3c, 1, sd), vc3 = apply(data3c, 1, var))

dataCext1 = data.frame(id = data_30mln_noX[[1]]$ensembl_gene_id,
                       mc = rowMeans(data1c), sdc = apply(data1c, 1, sd), vc = apply(data1c, 1, var),
                       who = list_of_libprepnames[[1]], what = " original")
dataCext2 = data.frame(id = data_30mln_noX[[1]]$ensembl_gene_id,
                       mc = rowMeans(data2c), sdc = apply(data2c, 1, sd), vc = apply(data2c, 1, var),
                       who = list_of_libprepnames[[2]], what = " original")
dataCext3 = data.frame(id = data_30mln_noX[[1]]$ensembl_gene_id,
                       mc = rowMeans(data3c), sdc = apply(data3c, 1, sd), vc = apply(data3c, 1, var),
                       who = list_of_libprepnames[[3]], what = " original")
dataCext1acc = dataCext1
dataCext1acc$vc = dataCext1acc$vc/(list_of_constsM[[1]])**2
dataCext1acc$what = "divided by QCC^2"
dataCext2acc = dataCext2
dataCext2acc$vc = dataCext2acc$vc/(list_of_constsM[[2]])**2
dataCext2acc$what = "divided by QCC^2"
dataCext3acc = dataCext3
dataCext3acc$vc = dataCext3acc$vc/(list_of_constsM[[3]])**2
dataCext3acc$what = "divided by QCC^2"

dataCext = rbind(dataCext1, dataCext2, dataCext3)
dataCextacc = rbind(dataCext1acc, dataCext2acc, dataCext3acc)

dataCexlALL = rbind(dataCext, dataCextacc)

dataCexlALLlogG1 = dataCexlALL[dataCexlALL$mc>=1, ]
dataCexlALLlogG1[, 2:4] = log(dataCexlALLlogG1[, 2:4], base=10)

dataCexlALLlogG10 = dataCexlALL[dataCexlALL$mc>=10, ]
dataCexlALLlogG10[, 2:4] = log(dataCexlALLlogG10[, 2:4], base=10)

#Fig.4C

libnreps = c(5,6,6)
reps_N_from30 = list(reps_6_from30[[1]][2:6], reps_6_from30[[2]], reps_6_from30[[3]])

sdPairs = lapply(1:3, function(l){
  combn(1:libnreps[[l]], 2, function(x){
    df = data_30mln_noX[[l]][, sort(c(1, reps_N_from30[[l]][x]*2, reps_N_from30[[l]][x]*2+1))]
    var_cov = sapply(1:nrow(df), function(i){var(c(df[i,2]+df[i,3], df[i,4]+df[i,5]))})
    cov = MeanCoverage(df)$meanCOV*2
    df_log = data.frame(cov_log=log(cov, base=10), var_cov_log=log(var_cov, base=10))
    df_log = df_log[!is.infinite(df_log$cov_log) & !is.infinite(df_log$var_cov_log),]
    sqrt_log_intercept = sqrt(10**lm(df_log, formula = (var_cov_log ~ offset(cov_log)))$coefficients)
    return(sqrt_log_intercept)
  })
})

slopesG1 = sapply(unique(dataCexlALLlogG1$what), function(what){
  sapply(unique(dataCexlALLlogG1$who), function(who){
    sqrt(10**(lm(dataCexlALLlogG1[dataCexlALLlogG1$who==who & dataCexlALLlogG1$what==what &
                                    dataCexlALLlogG1$vc!= -Inf, ], formula = (vc ~ offset(mc)))$coefficients))
  })
})
slopesG10 = sapply(unique(dataCexlALLlogG10$what), function(what){
  sapply(unique(dataCexlALLlogG10$who), function(who){
    sqrt(10**lm(dataCexlALLlogG10[dataCexlALLlogG10$who==who & dataCexlALLlogG10$what==what, ], formula = (vc ~ offset(mc)))$coefficients)
  })
})


CC_ov_sdp_df = data.frame(QCC = unlist(CC)[6:45],
                          sd_pair = unlist(sdPairs),
                          overdispersion_sqr = rep(slopesG10[,1], each=15)[6:45],
                          method = rep(unlist(list_of_libprepnames), each=15)[6:45])


#Fig.4A

df_2exp = lapply(1:3, function(i){
  data.frame(method = list_of_libprepnames[[i]],
             t(combn(reps_6_from30[[i]],2)),
             CC = CC[[i]])
})

df_1exp = lapply(1:3, function(i){
  do.call(rbind, lapply(1:ncol(combn(reps_6_from30[[i]],2)),
                        function(j){
                          pair1 = combn(reps_6_from30[[i]],2)[,j]
                          rest1 = reps_6_from30[[i]][! reps_6_from30[[i]] %in% pair1]
                          df = data.frame(X1=rep(pair1[1], each=ncol(combn(rest1, 2))),
                                          X2=rep(pair1[2], each=ncol(combn(rest1, 2))),
                                          X3=(combn(rest1, 2))[1,],
                                          X4=(combn(rest1, 2))[2,],
                                          CC_12 = CC[[i]][j])
                          df$CC_34 = sapply(1:nrow(df), function(k){
                            df_2exp[[i]][df_2exp[[i]]$X1==df[k, ]$X3 & df_2exp[[i]]$X2==df[k, ]$X4, ]$CC
                          })
                          df
                        }))
})

set.seed(3); df_1exp_15 = lapply(df_1exp, function(x){x[sample(nrow(x),15), ]})
set.seed(4); df_2exp_15 = lapply(1:ncol(combn(1:3,2)), function(c){
  m1 = combn(1:3,2)[1,c]
  m2 = combn(1:3,2)[2,c]
  data.frame(df_2exp[[m1]][sample(1:nrow(df_2exp[[m1]]),15), ],
             df_2exp[[m2]][sample(1:nrow(df_2exp[[m2]]),15), ])
})

interrep_exp_diff = lapply(1:3, function(df_1exp_n){
  lapply(1:15, function(i){
    data1 = df_1exp_15[[df_1exp_n]][i, ]
    df = data_30mln_noX[[df_1exp_n]][, c(1,
                                         sort(c(data1$X1*2, data1$X1*2+1, data1$X2*2, data1$X2*2+1)),
                                         sort(c(data1$X3*2, data1$X3*2+1, data1$X4*2, data1$X4*2+1)))]
    PerformBinTestAIAnalysisForTwoConditions_knownCC(df, vect1CondReps=1:2, vect2CondReps=3:4,
                                                     vect1CondRepsCombsCC=data1$CC_12, vect2CondRepsCombsCC=data1$CC_34,
                                                     Q=0.95, thr=10, thrUP=NA, thrType="each", minDifference=NA)
  })
})

betweenexp_diff = lapply(1:3, function(df_2exp_n){
  m1 = combn(1:3, 2)[1, df_2exp_n]
  m2 = combn(1:3, 2)[2, df_2exp_n]
  lapply(1:15, function(i){
    data2 = df_2exp_15[[df_2exp_n]][i, ]
    df1 = data_30mln_noX[[m1]][, c(1, sort(c(data2$X1*2, data2$X1*2+1, data2$X2*2, data2$X2*2+1)))]
    df2 = data_30mln_noX[[m2]][, c(1, sort(c(data2$X1.1*2, data2$X1.1*2+1, data2$X2.1*2, data2$X2.1*2+1)))]
    df = merge(df1, df2, by="ensembl_gene_id")
    head(df)
    PerformBinTestAIAnalysisForTwoConditions_knownCC(df, vect1CondReps=1:2, vect2CondReps=3:4,
                                                     vect1CondRepsCombsCC=data2$CC, vect2CondRepsCombsCC=data2$CC.1,
                                                     Q=0.95, thr=10, thrUP=NA, thrType="each", minDifference=NA)
  })
})

interrep_exp_diff_list = lapply(interrep_exp_diff, function(m_list){
  do.call(rbind, lapply(1:length(m_list), function(x){
    df = na.omit(m_list[[x]][, c("BT","BT_CC")])
    data.frame(
      BT_count = sum(df[, 1]),
      BTCC_count = sum(df[, 2]),
      count = nrow(df),
      BT_p = sum(df[, 1])/nrow(df),
      BTCC_p = sum(df[, 2])/nrow(df)
    )
  }))
})
interrep_exp_diff_df = do.call(rbind, lapply(1:3, function(i){
  rbind(
    data.frame(method = list_of_libprepnames[[i]],
               DiffASE_proportion = interrep_exp_diff_list[[i]]$BT_p,
               DiffASE_N = interrep_exp_diff_list[[i]]$BT_count,
               test = "BT",
               what = "1 Experiment"),
    data.frame(method = list_of_libprepnames[[i]],
               DiffASE_proportion = interrep_exp_diff_list[[i]]$BTCC_p,
               DiffASE_N = interrep_exp_diff_list[[i]]$BTCC_count,
               test = "BTCC",
               what = "1 Experiment")
  )
}))

betweenexp_diff_list = lapply(betweenexp_diff, function(m_list){
  do.call(rbind, lapply(1:length(m_list), function(x){
    df = na.omit(m_list[[x]][, c("BT","BT_CC")])
    data.frame(
      BT_count = sum(df[, 1]),
      BTCC_count = sum(df[, 2]),
      count = nrow(df),
      BT_p = sum(df[, 1])/nrow(df),
      BTCC_p = sum(df[, 2])/nrow(df)
    )
  }))
})

betweenexp_exp_diff_df = do.call(rbind, lapply(1:3, function(i){
  mm = combn(1:3, 2)[,i]
  rbind(
    data.frame(method = paste(list_of_libprepnames[[mm[1]]], "vs", list_of_libprepnames[[mm[2]]]),
               DiffASE_proportion = betweenexp_diff_list[[i]]$BT_p,
               DiffASE_N = betweenexp_diff_list[[i]]$BT_count,
               test = "BT",
               what = "2 Experiments"),
    data.frame(method = paste(list_of_libprepnames[[mm[1]]], "vs", list_of_libprepnames[[mm[2]]]),
               DiffASE_proportion = betweenexp_diff_list[[i]]$BTCC_p,
               DiffASE_N = betweenexp_diff_list[[i]]$BTCC_count,
               test = "BTCC",
               what = "2 Experiments")
  )
}))

df_fig1A <- rbind(betweenexp_exp_diff_df, interrep_exp_diff_df)
df_fig1A$method <- factor(df_fig1A$method,
                             levels = c("NEBNext (100ng) vs SMARTseq (10ng)",
                                        "NEBNext (100ng) vs SMARTseq (0.1ng)",
                                        "SMARTseq (10ng) vs SMARTseq (0.1ng)",
                                        "NEBNext (100ng)",
                                        "SMARTseq (10ng)",
                                        "SMARTseq (0.1ng)"))
df_fig1A$test <- factor(df_fig1A$test, labels = c("before correction", "QCC corrected"))
# Plotting


figure_4A <- ggplot(df_fig1A, aes(method, DiffASE_N, col= method)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~ test) +
  scale_color_manual(values=c("black","black","black", "royalblue1","maroon2","olivedrab3")) +
  ylab("Number of genes called differential") +
  xlab("") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "None", text = element_text(size=18))


df_fig4B <- rbind(res62_df_all,res61_df_all)
df_fig4B$method <- factor(df_fig4B$method, labels = c("1 replicate,\nno correction", "2 replicates,\nno correction", "2 replicates,\nQCC correction"))

figure_4B <- ggplot(df_fig4B, aes(x=experiment, y=FP_rate, col=experiment)) +
  geom_boxplot() +
  facet_grid(~ method) +
  xlab("") +
  ylab("False positive rate") +
  theme_bw() +
  scale_color_manual(values=c("royalblue1","olivedrab3","maroon2")) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), text = element_text(size=18), legend.position = "None")

cowplot::save_plot(figure_4A, file="~/Dropbox (Partners HealthCare)/replicates_ASE/manuscript/Figures/fig.4A_legend.pdf", base_height = 10, base_width = 16)


figure_4C <- grid.arrange(
  ggplot() +
    geom_boxplot(data=CC_ov_sdp_df, aes(overdispersion_sqr, QCC, col=method)) +
    geom_point(data=CC_ov_sdp_df, aes(overdispersion_sqr, QCC, col=method)) +
    theme_bw() +
    scale_color_manual(values=c("royalblue1","olivedrab3","maroon2")) +
    xlab("Sqrt of overdispersion") +
    theme(legend.position= "None", text = element_text(size=18)) ,
  ggplot() +
    geom_point(data=CC_ov_sdp_df, aes(sd_pair, QCC, col=method)) +
    theme_bw() +
    scale_color_manual(values=c("royalblue1","olivedrab3","maroon2")) +
    xlab("Counts sd") +
    theme(legend.position= "None", text = element_text(size=18)) ,
  ncol=2
)

figure_4D <- ggplot() +
  geom_point(data=dataCexlALLlogG1, aes(exp(mc), exp(vc), col = who), size=0.7, alpha=0.4) +
  geom_abline(slope = 1, linetype="dashed") +
  scale_y_log10() +
  scale_x_log10() +
  facet_grid(~what) +
  theme_bw() +
  stat_smooth(data=dataCexlALLlogG1, aes(exp(mc), exp(vc), group=who), col="black", method = "lm", formula = (y ~ offset(x)))+
  scale_color_manual(values=c("royalblue1","olivedrab3","maroon2")) +
  ylab("Dispersion of gene coverages") + xlab("Mean gene coverage") +
  theme(legend.position="None", text = element_text(size=18)) +
  coord_fixed()



PLT_fig4 = plot_grid(
  plot_grid(figure_4A, figure_4B, labels = c("A", "B"), rel_widths = c(0.8, 0.8)),
  plot_grid(figure_4C, figure_4D, labels = c("C", "D"), rel_widths = c(0.8, 0.8)),
  nrow = 2,
  scale = c(0.9, 0.9, 0.9, 0.9, 0.9)
)

cowplot::save_plot(PLT_fig4, file="~/Dropbox (Partners HealthCare)/replicates_ASE/manuscript/Figures/fig.4_v2.pdf", base_height = 10, base_width = 16)
