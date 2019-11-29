setwd("~/Dropbox (Partners HealthCare)/code/ASE/plot_figures")

library(tidyverse)
library(cowplot)
library(eulerr)
library(gridExtra)

source("../R/ASE_functions.R")
source("../R/PerformAIAnalysis_CC.v0.3.R")
source("../R/DownstreamASE.R")
source("../plot_figures/utilPlots.R")
source("../plot_figures/takedata30mln.R")

# reps_6_from30 = list(c(3,8,13,19,24,28),
#                      c(4,10,12,19,24,28),
#                      c(5,10,12,16,23,27))
reps_6_from30 = reps_6_fromXX
list_of_datas = data_XXmln_noX; list_of_datas[[1]] = data_XXmln_noX[[1]][, -c(2:11)]
list_of_consts = list(CC_noX[[1]], CC_noX[[2]], CC_noX[[3]])
list_of_constsM = list(mean(CC_noX[[1]][6:15]), mean(CC_noX[[2]]), mean(CC_noX[[3]]))
list_of_libprepnames = list("NEBNext (100ng)", "SMARTseq (10ng)", "SMARTseq (0.1ng)")
list_of_Nreps = list(6,6,6)

res62_df <-
  lapply(1:3, function(m){
    for_6_df <- CountsToAI(data_XXmln_noX[[m]], reps_6_from30[[m]], meth="mergedToProportion")$AI
    df <- sapply(1:length(list_of_consts[[m]]), function(j){
      pairs_sample = reps_6_from30[[m]][combn(list_of_Nreps[[m]], 2)[, j]]
      res62 <- PerformBinTestAIAnalysisForConditionNPointVect_knownCC(data_XXmln_noX[[m]],
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
    for_6_df <- CountsToAI(data_XXmln_noX[[m]], reps_6_from30[[m]], meth="mergedToProportion")$AI
    df <- sapply(1:list_of_Nreps[[m]], function(j){
      res61 <- PerformBinTestAIAnalysisForConditionNPointVect_knownCC(data_XXmln_noX[[m]],  #list_of_datas[[m]],
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


###-------------------------------------------------------------------------------------------------
###--------------------------------FIG_4D-----------------------------------------------------------
###-------------------------------------------------------------------------------------------------

data1c = sapply(reps_6_from30[[1]][2:6], function(i){rowSums(data_XXmln_noX[[1]][, c(i*2,i*2+1)])})
data2c = sapply(reps_6_from30[[2]][1:6], function(i){rowSums(data_XXmln_noX[[2]][, c(i*2,i*2+1)])})
data3c = sapply(reps_6_from30[[3]][1:6], function(i){rowSums(data_XXmln_noX[[3]][, c(i*2,i*2)+1])})


dataC = data.frame(id = data_XXmln_noX[[1]]$ID,
                   mc1 = rowMeans(data1c), sdc1 = apply(data1c, 1, sd), vc1 = apply(data1c, 1, var),
                   mc2 = rowMeans(data2c), sdc2 = apply(data2c, 1, sd), vc2 = apply(data2c, 1, var),
                   mc3 = rowMeans(data3c), sdc3 = apply(data3c, 1, sd), vc3 = apply(data3c, 1, var))

dataCext1 = data.frame(id = data_XXmln_noX[[1]]$ID,
                       mc = rowMeans(data1c), sdc = apply(data1c, 1, sd), vc = apply(data1c, 1, var),
                       who = list_of_libprepnames[[1]], what = " original")
dataCext2 = data.frame(id = data_XXmln_noX[[1]]$ID,
                       mc = rowMeans(data2c), sdc = apply(data2c, 1, sd), vc = apply(data2c, 1, var),
                       who = list_of_libprepnames[[2]], what = " original")
dataCext3 = data.frame(id = data_XXmln_noX[[1]]$ID,
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

###-------------------------------------------------------------------------------------------------
###--------------------------------FIG_4E-----------------------------------------------------------
###-------------------------------------------------------------------------------------------------

libnreps = c(5,6,6)
reps_N_from30 = list(reps_6_from30[[1]][2:6], reps_6_from30[[2]], reps_6_from30[[3]])
#libnreps = c(6,6,6)
#reps_N_from30 = list(reps_6_from30[[1]], reps_6_from30[[2]], reps_6_from30[[3]])

sdPairs = lapply(1:3, function(l){
  combn(1:libnreps[[l]], 2, function(x){
    df = data_XXmln_noX[[l]][, sort(c(1, reps_N_from30[[l]][x]*2, reps_N_from30[[l]][x]*2+1))]
    var_cov = sapply(1:nrow(df), function(i){var(c(df[i,2]+df[i,3], df[i,4]+df[i,5]))})
    cov = MeanCoverage(df)$meanCOV*2
    df_log = data.frame(cov_log=log(cov, base=10), var_cov_log=log(var_cov, base=10))
    df_log = df_log[!is.infinite(df_log$cov_log) & !is.infinite(df_log$var_cov_log),]
    sqrt_log_intercept = sqrt(10**lm(df_log, formula = (var_cov_log ~ offset(cov_log)))$coefficients)
    return(sqrt_log_intercept)
  })
})

sdPairs_2 = lapply(1:3, function(l){
  combn(1:libnreps[[l]], 2, function(x){
    df = data_XXmln_noX[[l]][, sort(c(1, reps_N_from30[[l]][x]*2, reps_N_from30[[l]][x]*2+1))]
    df = na.omit(ThresholdingCounts(df, thr=10))
    var_cov = sapply(1:nrow(df), function(i){var(c(df[i,2]+df[i,3], df[i,4]+df[i,5]))})
    cov = MeanCoverage(df)$meanCOV*2
    var_mean_ratio = exp(mean(log(var_cov/cov)[!is.infinite(log(var_cov/cov))]))
    return(var_mean_ratio)
  })
})

slopes_generator = function(df_log10){
  sapply(unique(df_log10$what), function(what){
    sapply(unique(df_log10$who), function(who){
      sqrt(10**(lm(df_log10[df_log10$who==who & df_log10$what==what & df_log10$mc!= -Inf & df_log10$vc!= -Inf, ],
                   formula = (vc ~ offset(mc)))$coefficients))
    })
  })
}

slopesG1  = slopes_generator(dataCexlALLlogG1)
slopesG10 = slopes_generator(dataCexlALLlogG10)

var_mean_ratio = sapply(unlist(list_of_libprepnames), function(who){
  df = dataCext[dataCext$mc >= 10 & dataCext$who == who,]
  exp(mean(log(df$vc / df$mc)[!is.infinite(log(df$vc / df$mc))]))
})


CC_ov_sdp_df = data.frame(QCC = unlist(CC)[6:45],
                          sd_pair = unlist(sdPairs),
                          overdispersion_sqr = rep(slopesG10[,1], each=15)[6:45],
                          method = rep(unlist(list_of_libprepnames), each=15)[6:45])
CC_ov_sdp_df_2 = data.frame(QCC = unlist(CC)[6:45],
                            sd_pair = unlist(sdPairs_2),
                            overdispersion_sqr = rep(sqrt(var_mean_ratio), each=15)[6:45],
                            method = rep(unlist(list_of_libprepnames), each=15)[6:45])
# CC_ov_sdp_df = data.frame(QCC = unlist(CC),
#                           sd_pair = unlist(sdPairs),
#                           overdispersion_sqr = rep(slopesG10[,1], each=15),
#                           method = rep(unlist(list_of_libprepnames), each=15))

###-------------------------------------------------------------------------------------------------
###--------------------------------FIG_4A-----------------------------------------------------------
###-------------------------------------------------------------------------------------------------

df_2exp = lapply(1:3, function(i){
  data.frame(method = list_of_libprepnames[[i]],
             t(combn(reps_6_from30[[i]],2)),
             CC = CC_noX[[i]])
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
                                          CC_12 = CC_noX[[i]][j])
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
    df = data_XXmln_noX[[df_1exp_n]][, c(1,
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
    df1 = data_XXmln_noX[[m1]][, c(1, sort(c(data2$X1*2, data2$X1*2+1, data2$X2*2, data2$X2*2+1)))]
    df2 = data_XXmln_noX[[m2]][, c(1, sort(c(data2$X1.1*2, data2$X1.1*2+1, data2$X2.1*2, data2$X2.1*2+1)))]
    df = merge(df1, df2, by="ID")
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

df_fig1A <- rbind(interrep_exp_diff_df, betweenexp_exp_diff_df)
df_fig1A$method <- factor(df_fig1A$method,
                             levels = c("NEBNext (100ng)",
                                        "SMARTseq (10ng)",
                                        "SMARTseq (0.1ng)",
                                        "NEBNext (100ng) vs SMARTseq (10ng)",
                                        "NEBNext (100ng) vs SMARTseq (0.1ng)",
                                        "SMARTseq (10ng) vs SMARTseq (0.1ng)"))
df_fig1A$test <- factor(df_fig1A$test, labels = c("before correction", "QCC corrected"))

###-------------------------------------------------------------------------------------------------
###-------------------------------------FIG_4C(b)---------------------------------------------------
###-------------------------------------------------------------------------------------------------

CalculatePairConcordance <- function(data, cc, expname, nrep=6){

  pdis_CC =
    sapply(1:ncol(combn(nrep, 2)), function(j){
      #sample2reps30mln =(combn(6, 2)[, j]-1)*5 + sample(1:5, 2, replace=T)
      sample2reps30mln = (combn(nrep, 2)[, j])
      CC_sample2_30mln = cc
      data_XXmln_sample2 = data[, sort(c(1, sample2reps30mln*2, sample2reps30mln*2+1))]

      if(cc == 1) {
        p = GetPercDiffForExperiment_BT(df=data_XXmln_sample2, CC_df=CC_sample2_30mln, exp_name="any")
      } else {
        p = GetPercDiffForExperiment_BT_CC(df=data_XXmln_sample2, CC_df=CC_sample2_30mln, exp_name="any")
      }

      return(p)
    })

  df = data.frame(Pc=as.vector(pdis_CC[1,]), libprep=expname, what="corrected", t(combn(nrep, 2)))

  return(df)
}

Calculate62Desagreement <- function(data, cc, expname, nrep=6){

  df <- sapply(1:ncol(combn(nrep, 2)), function(j){
    #pairs_sample = (combn(nrep, 2)[, j]-1)*5 + sample(1:5, 2, replace=T)
    pairs_sample = (combn(nrep, 2)[, j])
    for_6_df <- CountsToAI(data, meth="mergedToProportion")$AI
    res62 <- PerformBinTestAIAnalysisForConditionNPointVect_knownCC(data, vectReps=pairs_sample,
                                                                    vectRepsCombsCC = cc,
                                                                    ptVect = for_6_df,
                                                                    thr=10)
    FP_BTCC <- sum(res62$BT_CC, na.rm = T)
    FP_BT <- sum(res62$BT, na.rm = T)
    all_BTCC <- length(res62$BT_CC) - sum(is.na(res62$BT_CC))
    all_BT <- length(res62$BT) - sum(is.na(res62$BT))
    return(c(FP_BTCC, FP_BT, all_BTCC, all_BT))
  })
  df <- data.frame(t(df), t(combn(nrep, 2)))
  colnames(df) <- c("FP_BTCC", "FP_BT", "all_BTCC", "all_BT", "X1", "X2")
  df$FP_BTCC_rate <- df$FP_BTCC / df$all_BTCC
  df$FP_BT_rate <- df$FP_BT / df$all_BT
  return(df)

  df$experiment <- expname
  res62_df <- df[,c(9,7:8)] %>% gather(key="method", value = "FP_rate", -experiment)

  return(res62_df)
}


CC_sim = lapply(seq(1,3.5,0.1), function(x){c(x)})

NEB_data_XXmln_noX = data_XXmln_noX[[1]][, sort(c(1, reps_6_from30[[1]]*2, reps_6_from30[[1]]*2+1))]
NEB_data_XXmln_noX_no1 = NEB_data_XXmln_noX[, -(2:3)]
SM10ng_data_XXmln_noX = data_XXmln_noX[[2]][, sort(c(1, reps_6_from30[[2]]*2, reps_6_from30[[2]]*2+1))]
SM100pg_data_XXmln_noX = data_XXmln_noX[[3]][, sort(c(1, reps_6_from30[[3]]*2, reps_6_from30[[3]]*2+1))]

data_list = list(NEB_data_XXmln_noX, NEB_data_XXmln_noX_no1, SM10ng_data_XXmln_noX, SM100pg_data_XXmln_noX)
nrep_list = c(6,5,6,6)
libprepname_list = c("NEBNext (100ng)", "NEBNext (100ng) 2-6 reps", "SMARTseq (10ng)", "SMARTseq (0.1ng)")

# res_62_CCsim_NEBno1 = lapply(CC_sim, function(sim_CC_i){
#   j = 2
#   data = data_list[[j]]
#   const = sim_CC_i
#   libprepname = libprepname_list[[j]]
#   res62 = Calculate62Desagreement(data, const, libprepname, nrep=nrep_list[[j]])
#   res62$CC = const
#   print(paste("CC =", sim_CC_i, "DONE"))
#   return(res62)
# })
#
# df62 = do.call(rbind, lapply(res_62_CCsim_NEBno1, function(x){x}))

res_62_CCsim_S100 = lapply(CC_sim, function(sim_CC_i){
  j = 4
  data = data_list[[j]]
  const = sim_CC_i
  libprepname = libprepname_list[[j]]
  res62 = Calculate62Desagreement(data, const, libprepname, nrep=nrep_list[[j]])
  res62$CC = const
  print(paste("CC =", sim_CC_i, "DONE"))
  return(res62)
})

df62_S100 = do.call(rbind, lapply(res_62_CCsim_S100, function(x){x}))


###-------------------------------------------------------------------------------------------------
###-------------------------------------FIG_4C(a)---------------------------------------------------
###-------------------------------------------------------------------------------------------------

# reppairs_NEBno1 = reps_6_from30[[1]][-1]
# df_1exp_CCsim_NEBno1 = lapply(CC_sim, function(sim_CC_i){
#   do.call(rbind, lapply(1:ncol(combn(reppairs_NEBno1,2)),
#           function(j){
#             pair1 = combn(reppairs_NEBno1,2)[,j]
#             rest1 = reppairs_NEBno1[! reppairs_NEBno1 %in% pair1]
#             df = data.frame(X1=rep(pair1[1], each=ncol(combn(rest1, 2))),
#                             X2=rep(pair1[2], each=ncol(combn(rest1, 2))),
#                             X3=(combn(rest1, 2))[1,],
#                             X4=(combn(rest1, 2))[2,],
#                             CC_12 = sim_CC_i,
#                             CC_34 = sim_CC_i)
#             df
#           }))
# })
#
# n_sample = 10
# set.seed(2); df_1exp_CCsim_NEBno1_15 = lapply(df_1exp_CCsim_NEBno1, function(x){x[sample(nrow(x),n_sample), ]})
#
# interrep_exp_CCsim_NEBno1_diff = lapply(df_1exp_CCsim_NEBno1_15, function(df){
#   lapply(1:n_sample, function(i){
#     data1 = df[i, ]
#     print(data1)
#     PerformBinTestAIAnalysisForTwoConditions_knownCC(NEB_data_XXmln_noX,
#                                                      vect1CondReps=((c(data1$X1,data1$X2)-1)%/%5+1),
#                                                      vect2CondReps=((c(data1$X3,data1$X4)-1)%/%5+1),
#                                                      vect1CondRepsCombsCC=data1$CC_12, vect2CondRepsCombsCC=data1$CC_34,
#                                                      Q=0.95, thr=10, thrUP=NA, thrType="each", minDifference=NA)
#   })
# })
#
# interrep_exp_CCsim_NEBno1_diff_list = lapply(interrep_exp_CCsim_NEBno1_diff, function(m_list){
#   do.call(rbind, lapply(1:length(m_list), function(x){
#     df = na.omit(m_list[[x]][, c("BT","BT_CC")])
#     data.frame(
#       BTCC_p = sum(df[, 2])/nrow(df),
#       BTCC_count = sum(df[, 2]),
#       count = nrow(df)
#     )
#   }))
# })
# interrep_CCsim_NEBno1_exp_diff_df = do.call(rbind, lapply(1:length(CC_sim), function(i){
#   data.frame(DiffASE_proportion = interrep_exp_CCsim_NEBno1_diff_list[[i]]$BTCC_p,
#              DiffASE_N = interrep_exp_CCsim_NEBno1_diff_list[[i]]$BTCC_count,
#              QCC = CC_sim[[i]])
# }))




reppairs_S100 = reps_6_from30[[3]]
df_1exp_CCsim_S100 = lapply(CC_sim, function(sim_CC_i){
  do.call(rbind, lapply(1:ncol(combn(reppairs_S100,2)),
                        function(j){
                          pair1 = combn(reppairs_S100,2)[,j]
                          rest1 = reppairs_S100[! reppairs_S100 %in% pair1]
                          df = data.frame(X1=rep(pair1[1], each=ncol(combn(rest1, 2))),
                                          X2=rep(pair1[2], each=ncol(combn(rest1, 2))),
                                          X3=(combn(rest1, 2))[1,],
                                          X4=(combn(rest1, 2))[2,],
                                          CC_12 = sim_CC_i,
                                          CC_34 = sim_CC_i)
                          df
                        }))
})

n_sample = 10
set.seed(2); df_1exp_CCsim_S100_15 = lapply(df_1exp_CCsim_S100, function(x){x[sample(nrow(x),n_sample), ]})

interrep_exp_CCsim_S100_diff = lapply(df_1exp_CCsim_S100_15, function(df){
  lapply(1:n_sample, function(i){
    data1 = df[i, ]
    print(data1)
    PerformBinTestAIAnalysisForTwoConditions_knownCC(SM100pg_data_XXmln_noX,
                                                     vect1CondReps=((c(data1$X1,data1$X2)-1)%/%5+1),
                                                     vect2CondReps=((c(data1$X3,data1$X4)-1)%/%5+1),
                                                     vect1CondRepsCombsCC=data1$CC_12, vect2CondRepsCombsCC=data1$CC_34,
                                                     Q=0.95, thr=10, thrUP=NA, thrType="each", minDifference=NA)
  })
})

interrep_exp_CCsim_S100_diff_list = lapply(interrep_exp_CCsim_S100_diff, function(m_list){
  do.call(rbind, lapply(1:length(m_list), function(x){
    df = na.omit(m_list[[x]][, c("BT","BT_CC")])
    data.frame(
      BTCC_p = sum(df[, 2])/nrow(df),
      BTCC_count = sum(df[, 2]),
      count = nrow(df)
    )
  }))
})
interrep_CCsim_S100_exp_diff_df = do.call(rbind, lapply(1:length(CC_sim), function(i){
  data.frame(DiffASE_proportion = interrep_exp_CCsim_S100_diff_list[[i]]$BTCC_p,
             DiffASE_N = interrep_exp_CCsim_S100_diff_list[[i]]$BTCC_count,
             QCC = CC_sim[[i]])
}))



###-------------------------------------------------------------------------------------------------
###-------------------------------------FIG_4C(c)---------------------------------------------------
###-------------------------------------------------------------------------------------------------

# res_22_CCsim_S100 = lapply(CC_sim, function(sim_CC_i){
#   j = 4
#   data = data_list[[j]]
#   const = sim_CC_i
#   libprepname = libprepname_list[[j]]
#   res22 = CalculatePairConcordance(data, const, libprepname, nrep=nrep_list[[j]])
#   res22$CC = const
#   print(paste("CC =", sim_CC_i, "DONE"))
#   return(res22)
# })
#
# df22_S100 = do.call(rbind, lapply(res_22_CCsim_S100, function(x){x}))
#
# figure_4Cc_SM100 = ggplot(df22_S100[df22_S100$what=="corrected",], aes(x=CC, y=Pc, group=CC)) +
#   geom_rect(xmin = min(CC_noX[[3]]), xmax = max(CC_noX[[3]]), ymin = -Inf, ymax = Inf,
#             fill=alpha("darkolivegreen1", 0.2)) +
#   geom_boxplot() +
#   geom_point() +
#   xlab("TEsted QCC") +
#   ylab("Concordance % \nbetween two replicates") +
#   theme_bw() +
#   theme(legend.position="None", text = element_text(size=18))

####
#### OFFTOP: simQCC for 17reps
####

# res_22_CCsim_S100_l1 = lapply(seq(1,7,0.1), function(sim_CC_i){
#   j = 4
#   chop = sample(1:6, 2)
#   data = data_list[[j]][, sort(c(1, chop*2, chop*2+1))]
#   const = sim_CC_i
#   libprepname = libprepname_list[[j]]
#   res22 = CalculatePairConcordance(data, const, libprepname, nrep=2)
#   res22$CC = const
#   print(paste("CC =", sim_CC_i, "DONE"))
#   print(res22)
#   return(res22)
# })
# df22_S100_l1 = do.call(rbind, lapply(res_22_CCsim_S100_l1, function(x){x}))
#
# res_22_CCsim_NEB_l1 = lapply(seq(1,7,0.1), function(sim_CC_i){
#   j = 2
#   chop = sample(1:5, 2)
#   data = data_list[[j]][, sort(c(1, chop*2, chop*2+1))]
#   const = sim_CC_i
#   libprepname = libprepname_list[[j]]
#   res22 = CalculatePairConcordance(data, const, libprepname, nrep=2)
#   res22$CC = const
#   print(paste("CC =", sim_CC_i, "DONE"))
#   print(res22)
#   return(res22)
# })
# df22_NEB_l1 = do.call(rbind, lapply(res_22_CCsim_NEB_l1, function(x){x}))
#
# res_22_CCsim_S10_l1 = lapply(seq(1,7,0.1), function(sim_CC_i){
#   j = 3
#   chop = sample(1:6, 2)
#   data = data_list[[j]][, sort(c(1, chop*2, chop*2+1))]
#   const = sim_CC_i
#   libprepname = libprepname_list[[j]]
#   res22 = CalculatePairConcordance(data, const, libprepname, nrep=2)
#   res22$CC = const
#   print(paste("CC =", sim_CC_i, "DONE"))
#   print(res22)
#   return(res22)
# })
# df22_S10_l1 = do.call(rbind, lapply(res_22_CCsim_S10_l1, function(x){x}))
#
#
# ggplot(rbind(df22_S100_l1[df22_S100_l1$what=="corrected",],
#              df22_NEB_l1[df22_NEB_l1$what=="corrected",],
#              df22_S10_l1[df22_S10_l1$what=="corrected",]),
#        aes(x=CC, y=Pc, group=libprep, color=libprep)) +
#   geom_rect(xmin = min(CC_noX[[1]][-c(1:5)]), xmax = max(CC_noX[[1]][-c(1:5)]), ymin = -Inf, ymax = Inf,
#             fill=alpha("skyblue1", 0.2), col="skyblue1") +
#   geom_rect(xmin = min(CC_noX[[3]]), xmax = max(CC_noX[[3]]), ymin = -Inf, ymax = Inf,
#             fill=alpha("darkolivegreen1", 0.2), col="darkolivegreen1") +
#   geom_rect(xmin = min(CC_noX[[2]]), xmax = max(CC_noX[[2]]), ymin = -Inf, ymax = Inf,
#             fill=alpha("lightpink", 0.2), col="lightpink") +
#   geom_point() +
#   scale_color_manual(values=c("royalblue1","olivedrab3","maroon2"))  +
#   xlab("TEsted QCC") +
#   ylab("Concordance % \nbetween two replicates") +
#   theme_bw() +
#   theme(legend.position="None", text = element_text(size=18))




###-------------------------------------------------------------------------------------------------
###===========================================PLOTS=================================================
###-------------------------------------------------------------------------------------------------


figure_4A <- ggplot(df_fig1A, aes(method, DiffASE_N, col= method)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(~ test) +
  scale_color_manual(values=c("royalblue1","maroon2","olivedrab3", "black","black","black")) +
  ylab("Number of genes called differential") +
  xlab("") +
  theme(axis.text.x=element_text(colour="white"), axis.ticks.x=element_blank(),
        legend.position = "None", text = element_text(size=18))


df_fig4B <- rbind(res62_df_all,res61_df_all)
df_fig4B$method <- factor(df_fig4B$method, labels = c("1 replicate,\nno correction", "2 replicates,\nno correction", "2 replicates,\nQCC correction"))
df_fig4B$experiment_n = sapply(df_fig4B$experiment, function(x){
  if(x=="NEBNext (100ng)") {
    1
  } else if(x=="SMARTseq (10ng)") {
    2
  } else {
    3
  }
})
df_fig4B$experiment = reorder(df_fig4B$experiment, df_fig4B$experiment_n)

figure_4B <- ggplot(df_fig4B, aes(x=experiment, y=FP_rate, col=experiment)) +
  geom_boxplot() +
  facet_grid(~ method) +
  xlab("") +
  ylab("False positive rate") +
  theme_bw() +
  scale_color_manual(values=c("royalblue1","maroon2","olivedrab3")) +
  theme(axis.text.x=element_text(colour="white"), axis.ticks.x=element_blank(),
        text = element_text(size=18), legend.position = "None")

cowplot::save_plot(figure_4A, file="fig.4A_legend.pdf", base_height = 10, base_width = 16)


figure_4E <- grid.arrange(
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
figure_4E_2 <- grid.arrange(
  ggplot() +
    geom_boxplot(data=CC_ov_sdp_df_2, aes(overdispersion_sqr, QCC, col=method)) +
    geom_point(data=CC_ov_sdp_df_2, aes(overdispersion_sqr, QCC, col=method)) +
    theme_bw() +
    scale_color_manual(values=c("royalblue1","olivedrab3","maroon2")) +
    xlab("Sqrt of overdispersion") +
    theme(legend.position= "None", text = element_text(size=18)) ,
  ggplot() +
    geom_point(data=CC_ov_sdp_df_2, aes(sd_pair, QCC, col=method)) +
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
  theme(legend.position="None", text = element_text(size=18))



# df_fig4Cb_NEB = df62
# figure_4Cb_NEB = ggplot(df_fig4Cb_NEB, aes(x=CC, y=FP_BTCC_rate, group=CC)) +
#   geom_rect(xmin = min(CC_noX[[1]][-c(1:5)]), xmax = max(CC_noX[[1]][-c(1:5)]), ymin = -Inf, ymax = Inf,
#             fill=alpha("skyblue1", 0.2)) +
#   geom_boxplot(color="grey", lwd=1.1) +
#   geom_point() +
#   xlab("Tested QCC") +
#   ylab("False positive rate") +
#   theme_bw() +
#   theme(legend.position="None", text = element_text(size=18))

# df_fig4Ca_NEB = interrep_CCsim_NEBno1_exp_diff_df
# figure_4Ca_NEB = ggplot(df_fig4Ca_NEB,
#                    aes(QCC, DiffASE_N, group=QCC)) +
#   geom_rect(xmin = min(CC_noX[[1]][-c(1:5)]), xmax = max(CC_noX[[1]][-c(1:5)]), ymin = -Inf, ymax = Inf,
#             fill=alpha("skyblue1", 0.2)) +
#   geom_boxplot(color="grey", lwd=1.1) +
#   geom_point() +
#   xlab("Tested QCC") +
#   ylab("Number of genes \ncalled differential") +
#   theme_bw() +
#   theme(legend.position="None", text = element_text(size=18))

df_fig4Cb_SM100 = df62_S100
figure_4Cb_SM100 = ggplot(df_fig4Cb_SM100, aes(x=CC, y=FP_BTCC_rate, group=CC)) +
  geom_rect(xmin = min(CC_noX[[3]]), xmax = max(CC_noX[[3]]), ymin = -Inf, ymax = Inf,
            fill=alpha("darkolivegreen1", 0.2)) +
  geom_boxplot(color="grey", lwd=1.1) +
  geom_point() +
  xlab("Tested QCC") +
  ylab("False positive rate") +
  theme_bw() +
  scale_y_log10() +
  theme(legend.position="None", text = element_text(size=18))

df_fig4Ca_SM100 = interrep_CCsim_S100_exp_diff_df
figure_4Ca_SM100 = ggplot(df_fig4Ca_SM100,
                   aes(QCC, DiffASE_N, group=QCC)) +
  geom_rect(xmin = min(CC_noX[[3]]), xmax = max(CC_noX[[3]]), ymin = -Inf, ymax = Inf,
            fill=alpha("darkolivegreen1", 0.2)) +
  geom_boxplot(color="grey", lwd=1.1) +
  geom_point() +
  xlab("Tested QCC") +
  ylab("Number of genes \ncalled differential") +
  theme_bw() +
  scale_y_log10() +
  theme(legend.position="None", text = element_text(size=18))

# log of zeroes handling:
df_fig4Cb_SM100_inf = df_fig4Cb_SM100; rplsr_b = min(df_fig4Cb_SM100_inf[df_fig4Cb_SM100_inf$FP_BTCC_rate!=0, "FP_BTCC_rate"]) * 0.6; df_fig4Cb_SM100_inf[df_fig4Cb_SM100_inf$FP_BTCC_rate==0, "FP_BTCC_rate"] = rplsr_b
figure_4Cb_SM100 = ggplot(df_fig4Cb_SM100_inf, aes(x=CC, y=FP_BTCC_rate, group=CC)) +
  geom_rect(xmin = min(CC_noX[[3]]), xmax = max(CC_noX[[3]]), ymin = -Inf, ymax = Inf,
            fill=alpha("darkolivegreen1", 0.2)) +
  geom_boxplot(color="grey", lwd=1.1, alpha = 0.3) +
  geom_point() +
  xlab("Tested QCC") +
  ylab("False positive rate") +
  theme_bw() +
  scale_y_log10(breaks=c(rplsr_b, 10**(-1:-4)), labels=c(0, formatC(10**(-1:-4),format="e",digits=0))) +
  theme(legend.position="None", text = element_text(size=18))

df_fig4Ca_SM100_inf = df_fig4Ca_SM100; rplsr = 10**(-0.2); df_fig4Ca_SM100_inf[df_fig4Ca_SM100_inf$DiffASE_N==0, "DiffASE_N"] = rplsr
figure_4Ca_SM100 = ggplot(df_fig4Ca_SM100_inf,
                          aes(QCC, DiffASE_N, group=QCC)) +
  geom_rect(xmin = min(CC_noX[[3]]), xmax = max(CC_noX[[3]]), ymin = -Inf, ymax = Inf,
            fill=alpha("darkolivegreen1", 0.2)) +
  geom_boxplot(color="grey", lwd=1.1,  alpha = 0.3) +
  geom_point() +
  xlab("Tested QCC") +
  ylab("Number of genes \ncalled differential") +
  theme_bw() +
  scale_y_log10(breaks=c(rplsr, 10**(0:3)), labels=c(0, 10**(0:3))) +
  theme(legend.position="None", text = element_text(size=18))
# -----------------------


figure_legends = ggplot(df_fig1A, aes(method, DiffASE_N,
                                      col= method)) +
  geom_point(size=3) +
  theme_bw() +
  scale_color_manual(values=c("black","black","black", "royalblue1","maroon2","olivedrab3")) +
  theme(axis.text.x=element_text(colour="white"), axis.ticks.x=element_blank(),
        legend.position = "top", text = element_text(size=18),
        legend.box.background = element_rect(colour = "grey"), legend.title = element_blank()) +
  guides(col = guide_legend(nrow = 3))

legend_plt = plot_grid(cowplot::get_legend(figure_legends), NULL, ncol=1, rel_heights = c(1, 0.3))



PLT_fig4_ABDE = plot_grid(
  plot_grid(figure_4A, NULL, figure_4B, labels = c("A", "", "B"), rel_widths = c(0.85, 0.1, 0.85), rel_heights = c(0.8, 0.8, 0.8), nrow=1),
  legend_plt,
  plot_grid(figure_4D, NULL, figure_4E, labels = c("D", "", "E"), rel_widths = c(0.85, 0.1, 0.85), rel_heights = c(0.8, 0.8, 0.8), nrow=1),
  nrow = 3, rel_heights = c(0.9, 0.15, 0.9),
  align = 'vh' #,scale = c(0.9, 0.9, 0.9, 0.9)
)

## plot_grid(
##   figure_4A, NULL, figure_4B, figure_4D, NULL, figure_4Ea,figure_4Eb,
##   rel_widths = c(0.85,0.1,0.85, 0.85,0.1,0.425,0.425),
##   rel_heights = c(0.9,0.9,0.9, 0.9,0.9,0.9,0.9),
##   labels = c("A", "", "B", "D", "", "E", ""),
##   nrow = 2, ncol = 3,
##   align = 'vh'
## )

# PLT_fig4_Cabc_SM100 = plot_grid(
#   figure_4Ca_SM100, NULL, figure_4Cb_SM100, NULL, figure_4Cc_SM100, labels = c("C", "", "", "", ""),
#   ncol = 1, align = "v",
#   rel_heights = c(0.8, 0.05, 0.8, 0.05, 0.8)
#   #scale = c(0.9, 0.9, 0.9)
# )
# PLT_fig4_ABDECabc_SM100 = plot_grid(
#   PLT_fig4_ABDE, NULL, PLT_fig4_Cabc_SM100,
#   ncol=3, align = 'vh',
#   rel_widths = c(1.47, 0.08, 0.65)
# )
#
# cowplot::save_plot(PLT_fig4_ABDECabc_SM100, file="fig.4_v4_SM100_ABDECabc.pdf", base_height = 11, base_width = 21)
# cowplot::save_plot(PLT_fig4_ABDECabc_SM100, file="fig.4_v4_SM100_ABDECabc.png", base_height = 11, base_width = 21)



PLT_fig4_Cab_SM100 = plot_grid(
  figure_4Ca_SM100, NULL, figure_4Cb_SM100, labels = c("C", "", ""),
  ncol = 1, align = "v",
  rel_heights = c(0.95, 0.05, 0.95)
  #scale = c(0.9, 0.9, 0.9)
)
PLT_fig4_ABDECab_SM100 = plot_grid(
  PLT_fig4_ABDE, NULL, PLT_fig4_Cab_SM100,
  ncol=3, align = 'h',
  rel_widths = c(1.47, 0.08, 0.65)
)

cowplot::save_plot(PLT_fig4_ABDECab_SM100, file="fig.4_v4_SM100_ABDECab.pdf", base_height = 11, base_width = 21)
cowplot::save_plot(PLT_fig4_ABDECab_SM100, file="fig.4_v4_SM100_ABDECab.png", base_height = 11, base_width = 21)


# PLT_fig4_Cab_NEB = plot_grid(
#   figure_4Ca_NEB, NULL, figure_4Cb_NEB, labels = c("E", "", "F"),
#   ncol = 1, align = "v",
#   rel_heights = c(0.9, 0.15, 0.9)
#   #scale = c(0.9, 0.9, 0.9)
# )
# PLT_fig4_ABDECab_NEB = plot_grid(
#   PLT_fig4_ABDE, NULL, PLT_fig4_Cab_NEB,
#   ncol=3, align = 'h',
#   rel_widths = c(1.6, 0.05, 0.5)
# )
#
# cowplot::save_plot(PLT_fig4_ABDECab_NEB, file="fig.4_v4_NEB_ABDECab.pdf", base_height = 11, base_width = 21)
# cowplot::save_plot(PLT_fig4_ABDECab_NEB, file="fig.4_v4_NEB_ABDECab.png", base_height = 11, base_width = 21)













# PLT_fig4_full_NEB = plot_grid(
#   plot_grid(figure_4A, figure_4B, figure_4E_NEB, labels = c("A", "B", "E"),
#             rel_widths = c(0.8, 0.8, 0.5), rel_heights = c(0.8, 0.8, 0.8),
#             ncol = 3, align="h"),
#   plot_grid(figure_4C, figure_4D, figure_4F_NEB, labels = c("C", "D", "F"),
#             rel_widths = c(0.9, 0.7, 0.5), rel_heights = c(1,1,1),
#             ncol = 3, align="h"),
#   nrow = 2,
#   scale = c(0.95, 0.94, 0.95, 0.95, 95, 95)
# )
# cowplot::save_plot(PLT_fig4_full_NEB, file="fig.4_v3_NEB.pdf", base_height = 11, base_width = 21)
# cowplot::save_plot(PLT_fig4_full_NEB, file="fig.4_v3_NEB.png", base_height = 11, base_width = 21)
#
# PLT_fig4_full_SM100 = plot_grid(
#   plot_grid(figure_4A, figure_4B, figure_4E_SM100, labels = c("A", "B", "E"),
#             rel_widths = c(0.8, 0.8, 0.5), rel_heights = c(0.8, 0.8, 0.8),
#             ncol = 3, align="h"),
#   plot_grid(figure_4C, figure_4D, figure_4F_SM100, labels = c("C", "D", "F"),
#             rel_widths = c(0.9, 0.7, 0.5), rel_heights = c(1,1,1),
#             ncol = 3, align="h"),
#   nrow = 2,
#   scale = c(0.95, 0.94, 0.95, 0.95, 95, 95)
# )
# cowplot::save_plot(PLT_fig4_full_SM100, file="fig.4_v3_SM100.pdf", base_height = 11, base_width = 21)
# cowplot::save_plot(PLT_fig4_full_SM100, file="fig.4_v3_SM100.png", base_height = 11, base_width = 21)




