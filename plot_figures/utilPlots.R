# ---------------------------------------------------------------------------------------
#                 UTIL FUNCTIONS FOR PLOTTING
# ---------------------------------------------------------------------------------------



#' Creates ...
#'
#' @param df_data
#' @param reppair
#' @param pairconst
#' @param libprepname
#' @return ????
#' @examples
#'
CreateForplotDF <- function(df_data, reppair, pairconst, libprepname) {
  df_1out = lapply(reppair, function(i){
    PerformBinTestAIAnalysisForConditionNPoint_knownCC(df_data, i, pairconst, thr=10)
  })
  df_1out

  df_bt = merge(df_1out[[1]][, c("ID", "BT", "BT_CC")], df_1out[[2]][, c("ID", "BT", "BT_CC")], by="ID")
  df_aicov = merge(merge(CountsToAI(df_data, reps=reppair[1],thr=10),
                         CountsToAI(df_data, reps=reppair[2], thr=10),
                         by="ensembl_gene_id"),
                   merge(MeanCoverage(df_data, reps=reppair[1], thr=10),
                         MeanCoverage(df_data, reps=reppair[2], thr=10),
                         by="ensembl_gene_id"),
                   by = "ensembl_gene_id")
  names(df_aicov)[1] = "ID"

  df_forbtplot = merge(df_bt, df_aicov, by = "ID")
  df_forbtplot$libprep = libprepname
  df_forbtplot
}

#' Creates ...
#'
#' @param forbtplots
#' @return ????
#' @examples
#'
CreateForplotDF_btNbtcc <- function(forbtplots) {
  df_forplot = do.call(rbind, forbtplots)

  df_forplot_bt = na.omit(data.frame(df_forplot[, !sapply(names(df_forplot), function(a){grepl("BT_CC.", a, fixed=TRUE)})],
                                     test = "binomial"))
  df_forplot_btcc = na.omit(data.frame(df_forplot[, !sapply(names(df_forplot), function(a){grepl("BT.", a, fixed=TRUE)})],
                                       test = "corrected binomial"))
  names(df_forplot_btcc) = names(df_forplot_bt)

  df_forplot_bt$BT.xyeq = (df_forplot_bt$BT.x & df_forplot_bt$BT.y | !df_forplot_bt$BT.x & !df_forplot_bt$BT.y)
  df_forplot_btcc$BT.xyeq = (df_forplot_btcc$BT.x & df_forplot_btcc$BT.y | !df_forplot_btcc$BT.x & !df_forplot_btcc$BT.y)

  list(BTBF = df_forplot_bt, BTBFCC = df_forplot_btcc)
}


#' Creates ...
#'
#' @param DF_forplot
#' @return ????
#' @examples
#
CreateForplotDF_btNbtcc_colorescapers <- function(DF_forplot){
  res = lapply(DF_forplot, function(df_forplot){
    percent_of_diff_color = do.call(rbind, lapply(unique(df_forplot$libprep), function(l){
      df = df_forplot[df_forplot$libprep==l,]
      data.frame(percentage = round(c(sum(df$BT.x & !df$BT.y)/sum(df$BT.x),
                                      sum(!df$BT.x & df$BT.y)/sum(!df$BT.x),
                                      sum(df$BT.y & !df$BT.x)/sum(df$BT.y),
                                      sum(!df$BT.y & df$BT.x)/sum(!df$BT.y),
                                      sum(df$BT.x & df$BT.y)/sum(df$BT.x | df$BT.y)
      )*100, 1),
      numberOfGenes = c(sum(df$BT.x), sum(!df$BT.x), sum(df$BT.y), sum(!df$BT.y), NA),
      who = c("y_color_not_like_x_division", "y_color_not_like_x_division", "x_color_not_like_y_division", "x_color_not_like_y_division", "concordance_rate"),
      BT.x = c(T, F, F, T, NA),
      BT.y = c(F, T, T, F, NA),
      libprep = l
      )
    }))
    percent_of_diff_color$P_color_escapers = paste(percent_of_diff_color$percentage, "%")
    percent_of_diff_color
  })
  names(res) = names(DF_forplot)
  res
}

#' Creates dataframe with all data needed for plotting (for Figure 1 - no correction)
#'
#' @param df Input dataframe
#' @param CC_df Dataframe with precalculated correction constants
#' @param exp_name Name of the experiment
#' @param exp_n Number of the experiment (order in the dataframe)
#' @return Dataframe for plotting
#' @examples
#
GetDataForExperiment_BT <- function(df, CC_df, exp_name, exp_n) {
  df <- CreateForplotDF_btNbtcc(list(CreateForplotDF(df[[exp_n]], 1:2,
                                                     pairconst=CC_df[[exp_n]],
                                                     exp_name)))
  df_BT = df$BTBF
  df_BT$BT.x = as.factor(df_BT$BT.x)
  levels(df_BT$BT.x) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")
  df_BT$BT.y = as.factor(df_BT$BT.y)
  levels(df_BT$BT.y) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")
  return(df_BT)
}

#' Creates dataframe with all data needed for plotting (for Figure 3 - with correction)
#'
#' @param df Input dataframe
#' @param CC_df Dataframe with precalculated correction constants
#' @param exp_name Name of the experiment
#' @param exp_n Number of the experiment (order in the dataframe)
#' @return Dataframe for plotting
#' @examples
#
GetDataForExperiment_BTCC <- function(df, CC_df, exp_name, exp_n) {
  df <- CreateForplotDF_btNbtcc(list(CreateForplotDF(df[[exp_n]], 1:2,
                                                     pairconst=CC_df[[exp_n]],
                                                     exp_name)))
  df_BT = df$BTBFCC
  df_BT$BT.x = as.factor(df_BT$BT.x)
  levels(df_BT$BT.x) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")
  df_BT$BT.y = as.factor(df_BT$BT.y)
  levels(df_BT$BT.y) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")
  return(df_BT)
}

#' Gets percents of disagreement for a given dataset
#'
#' @param df Input dataframe
#' @param CC_df Dataframe with precalculated correction constants
#' @param exp_name Name of the experiment
#' @param exp_n Number of the experiment (order in the dataframe)
#' @return Dataframe for plotting
#' @examples
#
GetPercDiffForExperiment_BT <- function(df, CC_df, exp_name) {
  df <- CreateForplotDF_btNbtcc(list(CreateForplotDF(df, 1:2,
                                                     pairconst=CC_df,
                                                     exp_name)))
  percent_of_diff_color_df = CreateForplotDF_btNbtcc_colorescapers(df)
  percent_of_diff_color_df$BTBF$test = "binomial"
  percent_of_diff_color_df_BT = percent_of_diff_color_df$BTBF

  percent_of_diff_color_df_BT$BT.x = as.factor(percent_of_diff_color_df_BT$BT.x)
  levels(percent_of_diff_color_df_BT$BT.x) = c("Balanced genes \n H0=0.5 Not Rejected", "Genes with AI \n H0=0.5 Rejected")

  #p = percent_of_diff_color_df_BT[percent_of_diff_color_df_BT$who == "y_color_not_like_x_division" & percent_of_diff_color_df_BT$BT.x == "Genes with AI \n H0=0.5 Rejected" |
  #                                  percent_of_diff_color_df_BT$who == "x_color_not_like_y_division" & percent_of_diff_color_df_BT$BT.x == "Balanced genes \n H0=0.5 Not Rejected", ]$percentage
  p = percent_of_diff_color_df_BT[percent_of_diff_color_df_BT$who == "concordance_rate", ]$percentage
  return(p)
}

#' Gets percents of disagreement for a given dataset
#'
#' @param df Input dataframe
#' @param CC_df Dataframe with precalculated correction constants
#' @param exp_name Name of the experiment
#' @param exp_n Number of the experiment (order in the dataframe)
#' @return Dataframe for plotting
#' @examples
#
GetPercDiffForExperiment_BT_CC <- function(df, CC_df, exp_name) {
  df <- CreateForplotDF_btNbtcc(list(CreateForplotDF(df, 1:2,
                                                     pairconst=CC_df,
                                                     exp_name)))
  percent_of_diff_color_df = CreateForplotDF_btNbtcc_colorescapers(df)
  percent_of_diff_color_df$BTBFCC$test = "corrected binomial"
  percent_of_diff_color_df_BTCC = percent_of_diff_color_df$BTBFCC

  percent_of_diff_color_df_BTCC$BT.x = as.factor(percent_of_diff_color_df_BTCC$BT.x)
  levels(percent_of_diff_color_df_BTCC$BT.x) = c("Balanced genes \n H0=0.5 Not Rejected", "Genes with AI \n H0=0.5 Rejected")

  #p = percent_of_diff_color_df_BTCC[percent_of_diff_color_df_BTCC$who == "y_color_not_like_x_division" & percent_of_diff_color_df_BTCC$BT.x == "Genes with AI \n H0=0.5 Rejected" |
  #                                    percent_of_diff_color_df_BTCC$who == "x_color_not_like_y_division" & percent_of_diff_color_df_BTCC$BT.x == "Balanced genes \n H0=0.5 Not Rejected", ]$percentage
  p = percent_of_diff_color_df_BTCC[percent_of_diff_color_df_BTCC$who == "concordance_rate", ]$percentage
  return(p)
}






