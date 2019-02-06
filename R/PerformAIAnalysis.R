#
# PERFORM DIFF AI ANALYSIS ON 2 CONDITIONS
# _______________________________________________________________________________________
# TODO
#'
#'* add correction of CI on mean(AI_gene)
#'* add correction on overdispersion
#'* add bins->limits
#'* more smart lm needed
#'* think about renaming columns to rep1_ref and so on, it's bad

# _______________________________________________________________________________________

options(stringsAsFactors = FALSE)
# _______________________________________________________________________________________

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PERFORM DIFF CI(AI) ANALYSIS
# ---------------------------------------------------------------------------------------

PerformDAIQuantilesAnalysis <- function(inDF, vectReps, condName="Condition", 
                                       Q=0.95, EPS=1.3, thr=NA){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param condName An optional parameter; one-word name for condition
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param fullOUT Set true if you want full output with all computationally-internal dfs
  #' @return A table of quantiles in coverage bins 
  #' @examples
  #'
  
  # Take subtable:
  dfCondition <- inDF[, sort(c(1, vectReps*2, vectReps*2+1))]
  
  # Create pairvise AI differences for all techreps pairs:
  deltaAIPairwiseDF <- CreateMergedDeltaAIPairwiseDF(dfCondition, what=condName, thr=thr)
  deltaAIPairwiseDF$group <- paste(condName, deltaAIPairwiseDF$ij)
  
  # Count quantiles for Mean Coverage bins:
  observedQuantilesDF <- do.call(rbind,
                                 lapply(unique(deltaAIPairwiseDF$group),
                                        function(gr){
                                          df  <- deltaAIPairwiseDF[deltaAIPairwiseDF$group == gr, ]
                                          res <- CreateObservedQuantilesDF(df,
                                                                           P=Q, ep=EPS, logbase=T,
                                                                           coverageLimit=quantile(deltaAIPairwiseDF$MeanCov, 0.995),
                                                                           group=gr)
                                        }
                                 )
  )
  observedQuantilesDF$condition <- condName
  observedQuantilesDF$ij <- sapply(as.character(observedQuantilesDF$group),
                                   function(x){paste(unlist(strsplit(x, ' '))[2:4], collapse=' ')})
  
  return(observedQuantilesDF)
}

PerformCIAIAnalysis <- function(inDF, vectReps, condName="Condition", 
                                Q=0.95, EPS=1.3, thr=NA, fullOUT=F){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param condName An optional parameter; one-word name for condition
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param fullOUT Set true if you want full output with all computationally-internal dfs
  #' @return A table of gene names, AIs + CIs
  #' @examples
  #'
  
  # Take subtable:
  dfCondition <- inDF[, sort(c(1, vectReps*2, vectReps*2+1))]
  
  observedQuantilesDF <- PerformDAIQuantilesAnalysis(inDF, vectReps, condName, Q, EPS, thr)
  
  # Count intercepts:
  linIntercepts <- data.frame(condition = condName, 
                              ij = unique(observedQuantilesDF$ij),
                              linInt = as.double(sapply(unique(observedQuantilesDF$ij), function(x){ 
                                FitLmIntercept(observedQuantilesDF[observedQuantilesDF$ij == x, ], 
                                               binNObs=30, morethan=10, logoutput=F)
                              })))
                          
  # -------------------------------------------------------------------------------------
  # Concerns: Meancoverage is not accurate measure.
  #           Overdispersion itself should be enough, in case of 5aza
  #           it's even worse because the different level of coverage.
  # -------------------------------------------------------------------------------------
  
  # Calculate AI CIs:
  
  QCI <- CreateCIforAI(linIntercepts, dfCondition, thr=thr)
  
  if (!fullOUT){
    return(QCI)
  } else {
    return(list(AICI = QCI,
                FULL_OUT = list(observedQuantiles = observedQuantilesDF,
                                intercepts = linIntercepts)))
  }
}


# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PERFORM DIFF AI ANALYSIS ON 2 CONDITIONS
# ---------------------------------------------------------------------------------------
PerformDiffAIAnalysisFor2Conditions <- function(inDF, vect1CondReps, vect2CondReps,
                                                    cond1Name="Condition1", cond2Name="Condition2",
                                                    Q=0.95, EPS=1.3, thr=NA, fullOUT=F){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for each condition
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vect1CondReps A vector (>=2) of replicate numbers that should be considered as first condition's tech reps
  #' @param vect2CondReps A vector (>=2) of replicate numbers that should be considered as second condition's tech reps
  #' @param cond1Name An optional parameter; one-word name for condition 1
  #' @param cond2Name An optional parameter; one-word name for condition 2
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param fullOUT Set true if you want full output with all computationally-internal dfs.
  #' @return A table of gene names, AIs + CIs for each condition, classification into genes demonstrating differential AI and those that don't
  #' @examples
  #'
  
  # Bonferroni correction:
  Q <- 1 - (1-Q)/nrow(inDF)
  
  if (!fullOUT){
    QCI <- data.frame(inDF[, 1], 
                      cond1Name, 
                      PerformCIAIAnalysis(inDF, vect1CondReps, cond1Name, Q, EPS, thr, fullOUT=F)[, -1], 
                      cond2Name, 
                      PerformCIAIAnalysis(inDF, vect2CondReps, cond2Name, Q, EPS, thr, fullOUT=F)[, -1]) 
  } else {
    OUT <- list(cond1 = PerformCIAIAnalysis(inDF, vect1CondReps, cond1Name, Q, EPS, thr, fullOUT=T),
                cond2 = PerformCIAIAnalysis(inDF, vect2CondReps, cond2Name, Q, EPS, thr, fullOUT=T))
    QCI <- data.frame(inDF[, 1], 
                      cond1Name, OUT$cond1$AICI[, -1], 
                      cond2Name, OUT$cond2$AICI[, -1]) 
  }
  names(QCI) = c("ID",
                 "condition1", 
                 "meanCov1", "meanAI1", "pm1", "meanAI1Low", "meanAI1High",
                 "condition2",
                 "meanCov2", "meanAI2", "pm2", "meanAI2Low", "meanAI2High")
  
  # Find intersecting intervals > call them FALSE
  
  QCI$diffAI <- !(QCI$meanAI1Low < QCI$meanAI2Low & QCI$meanAI1High >= QCI$meanAI2Low |
                    QCI$meanAI1Low >= QCI$meanAI2Low & QCI$meanAI1Low <= QCI$meanAI2High)
  
  if (!fullOUT){
    return(QCI)
  } else {
    return(list(AICI = QCI,
                FULL_OUT = OUT))
  }
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PERFORM DIFF AI ANALYSIS ON CONDITION AND POINT ESTIMATION
# ---------------------------------------------------------------------------------------

PerformDiffAIAnalysisForConditionNPoint <- function(inDF, vectReps, condName="Condition", pt = 0.5,
                                                Q=0.95, EPS=1.3, thr=NA, fullOUT=F){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for condition + point estimate to compare
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param condName An optional parameter; one-word name for condition
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param fullOUT Set true if you want full output with all computationally-internal dfs
  #' @return A table of gene names, AIs + CIs, classification into genes demonstrating differential from point estimate AI and those that don't
  #' @examples
  #'
  
  # Bonferroni correction:
    Q <- 1 - (1-Q)/nrow(inDF)
  
  if (!fullOUT){
    QCI <- PerformCIAIAnalysis(inDF, vectReps, condName, Q, EPS, thr, fullOUT=F)
  } else {
    OUT <- PerformCIAIAnalysis(inDF, vectReps, condName, Q, EPS, thr, fullOUT=T)
    QCI <- OUT$AICI
  }
  
  # Find intersecting intervals > call them FALSE
  
  QCI$diffAI <- !(QCI$meanAILow <= pt & QCI$meanAIHigh >= pt)
  
  if (!fullOUT){
    return(QCI)
  } else {
    return(list(AICI = QCI,
                FULL_OUT = OUT))
  }
}



