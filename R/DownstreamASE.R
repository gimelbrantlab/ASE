# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PROCESS TABLES WITH COUNTs/AIs, FIND MAE GENES
# ---------------------------------------------------------------------------------------

DiffAIplusFig <- function(inDF, vect1CondReps_expA, vect1CondReps_expB,
                          vect1CondReps_expA_name, vect1CondReps_expB_name,
                          thr_coverage=NA, minDifference=0.1, lm=T) {
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vect1CondReps_expA A vector (>=2) of replicate numbers that should be considered as first condition's tech reps
  #' @param vect1CondReps_expB A vector (>=2) of replicate numbers that should be considered as second condition's tech reps
  #' @param vect1CondReps_expA_name Name of the first condition
  #' @param vect1CondReps_expB_name Name of the second condition
  #' @param thr_coverage optional parameter, specifies coverage threshold
  #' @param minDifference minimal difference between AIs in 2 conditions to call them differentially allelically expressed
  #' @return inDF table with 2 additional columns classifying genes into differentially allelically expressed (DAE) and not + figure
  #' @examples
  #'
  df_tmp <- PerformDiffAIAnalysisFor2Conditions(inDF,
                                                vect1CondReps = vect1CondReps_expA,
                                                vect2CondReps = vect1CondReps_expB,
                                                Q=0.95,
                                                thr=thr_coverage,
                                                EPS=1.3,
                                                fullOUT=F,
                                                minDifference=minDifference)

  fig_compare_tmp <- ggplot(df_tmp, aes(x = meanAI_1, y = meanAI_2, col = DAE)) +
    geom_point(size=0.5) +
    theme_bw() +
    xlab(paste0("AI, ",vect1CondReps_expA_name)) +
    ylab(paste0("AI, ",vect1CondReps_expB_name)) +
    scale_color_manual(name="Differential AI", labels=c("FALSE", "TRUE"), values=c("gray", "red")) +
    coord_fixed() +
    theme(legend.position = "None")
  if (lm)
    fig_compare_tmp <- fig_compare_tmp + geom_smooth(method="lm")

  return(list(df_tmp, fig_compare_tmp))
}



findMAE <- function(x) {
  #'
  #' @param x Vector with ASE classification per clone
  #' @return overall MAE status among clones (by Gendrel classification)
  #' @examples
  #'
  mon <- 0
  if (length(x)==1) {
    mon=NA
  }
  else {
    not_nm_count <- length(x) - sum(x=="nd")
    if (not_nm_count==0) mon="nd"
    else if ((sum(x=="CAST_monoallelic")+sum(x=="CAST_biased")==not_nm_count)|(sum(x=="129_monoallelic")+sum(x=="129_biased")==not_nm_count)) mon="gen_sk"
    else if ((sum(x=="CAST_monoallelic")>0)|(sum(x=="129_monoallelic")>0)) mon="monoallelic"
    else if ((sum(x=="CAST_biased")>0)|(sum(x=="129_biased")>0)) mon="biased"
    else if (sum(x=="biallelic")==not_nm_count) mon="biallelic"
    else mon="other"
  }
  return(mon)
}

isMAE_test_CI <- function(x) {
  #'
  #' @param x Vector two entries: first is the output of CIs test (T if genes are differentially monoallelically expressed), second is absolute AI value
  #' @return overall MAE status among clones (by Gendrel classification)
  #' @examples
  #'
  thr <- x[3]
  if ((is.na(x[1]))|(is.na(x[2]))) {
    return("nd")
  }
  else {
    if ((x[1])&((x[2]>=thr))) {
      return("129_monoallelic")
    }
    else if ((x[1])&((x[2]<=(1-thr)))) {
      return("CAST_monoallelic")
    }
    else if ((x[1])&((x[2]>=0.5))) {
      return("129_biased")
    }
    else if ((x[1])&((x[2]<0.5))) {
      return("CAST_biased")
    }
    else {return("biallelic")}
  }
}
