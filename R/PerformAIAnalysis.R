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
                                        Q=0.95, EPS=1.3, thr=NA, thrUP=NA, thrType="each"){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param condName An optional parameter; one-word name for condition
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @param fullOUT Set true if you want full output with all computationally-internal dfs
  #' @return A table of quantiles in coverage bins
  #' @examples
  #'

  # Take subtable:
  dfCondition <- inDF[, sort(c(1, vectReps*2, vectReps*2+1))]

  # Create pairvise AI differences for all techreps pairs:
  deltaAIPairwiseDF <- CreateMergedDeltaAIPairwiseDF(df=dfCondition, what=condName, thr=thr, thrUP=thrUP, thrType=thrType)
  deltaAIPairwiseDF$group <- paste(condName, deltaAIPairwiseDF$ij)

  # Count quantiles for Mean Coverage bins:
  observedQuantilesDF <- do.call(rbind,
                                 lapply(unique(deltaAIPairwiseDF$group),
                                        function(gr){
                                          df  <- deltaAIPairwiseDF[deltaAIPairwiseDF$group == gr, ]
                                          res <- CreateObservedQuantilesDF(df=df,
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
                                Q=0.95, EPS=1.3, thr=NA, thrUP=NA, thrType="each",
                                fullOUT=F){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param condName An optional parameter; one-word name for condition
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @param fullOUT Set true if you want full output with all computationally-internal dfs
  #' @return A table of gene names, AIs + CIs
  #' @examples
  #'

  # Take subtable:
  dfCondition <- inDF[, sort(c(1, vectReps*2, vectReps*2+1))]

  observedQuantilesDF <- PerformDAIQuantilesAnalysis(inDF=inDF, vectReps=vectReps, condName=condName,
                                                     Q=Q, EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType)

  # Count intercepts:
  if(is.na(thr)) {
    morethen = 10
  } else {
    morethen = thr
  }
  linIntercepts <- data.frame(condition = condName,
                              ij = unique(observedQuantilesDF$ij),
                              linInt = as.double(sapply(unique(observedQuantilesDF$ij), function(x){
                                FitLmIntercept(inDF=observedQuantilesDF[observedQuantilesDF$ij == x, ],
                                               binNObs=30, morethan=morethen, logoutput=F)
                              })))

  # -------------------------------------------------------------------------------------
  # Concerns: Meancoverage is not accurate measure.
  #           Overdispersion itself should be enough, in case of 5aza
  #           it's even worse because the different level of coverage.
  # -------------------------------------------------------------------------------------

  # Calculate AI CIs:

  QCI <- CreateCIforAI(dfInt=linIntercepts, dfCounts=dfCondition,
                       condName=condName, thr=thr, thrUP, thrType)

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
                                                Q=0.95, EPS=1.3, thr=NA, thrUP=NA, thrType="each",
                                                minDifference=0.1,
                                                fullOUT=F){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for each condition
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vect1CondReps A vector (>=2) of replicate numbers that should be considered as first condition's tech reps
  #' @param vect2CondReps A vector (>=2) of replicate numbers that should be considered as second condition's tech reps
  #' @param cond1Name An optional parameter; one-word name for condition 1
  #' @param cond2Name An optional parameter; one-word name for condition 2
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @param fullOUT Set true if you want full output with all computationally-internal dfs.
  #' @return A table of gene names, AIs + CIs for each condition, classification into genes demonstrating differential AI and those that don't
  #' @examples
  #'

  # Bonferroni correction:
  Q <- 1 - (1-Q)/nrow(inDF)

  if (!fullOUT){
    QCI <- data.frame(inDF[, 1],
                      PerformCIAIAnalysis(inDF=inDF, vectReps=vect1CondReps, condName=cond1Name,
                                          Q=Q, EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType, fullOUT=F)[, -1],
                      PerformCIAIAnalysis(inDF=inDF, vectReps=vect2CondReps, condName=cond2Name,
                                          Q=Q, EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType, fullOUT=F)[, -1])
  } else {
    OUT <- list(cond1 = PerformCIAIAnalysis(inDF=inDF, vectReps=vect1CondReps, condName=cond1Name,
                                            Q=Q, EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType, fullOUT=T),
                cond2 = PerformCIAIAnalysis(inDF=inDF, vectReps=vect2CondReps, condName=cond2Name,
                                            Q=Q, EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType, fullOUT=T))
    QCI <- data.frame(inDF[, 1],
                      OUT$cond1$AICI[, -1],
                      OUT$cond2$AICI[, -1])
  }
  names(QCI) = c("ID", paste0(c("condition", "meanCov", "meanAI",
                                "pm", "meanAILow", "meanAIHigh"),
                              '_', rep(1:2, each=6))
                )

  # Find intersecting intervals > call them FALSE (non-rejected H_0)
  QCI$diffAI <- !(QCI$meanAILow_1 < QCI$meanAILow_2 & QCI$meanAIHigh_1 >= QCI$meanAILow_2 |
                  QCI$meanAILow_1 >= QCI$meanAILow_2 & QCI$meanAILow_1 <= QCI$meanAIHigh_2)
  QCI$DAE <- (QCI$diffAI & (abs(QCI$meanAI_1 - QCI$meanAI_2) >= minDifference))

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
                                                    Q=0.95, EPS=1.3, thr=NA, thrUP=NA, thrType="each",
                                                    minDifference=0.1,
                                                    fullOUT=F){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for condition + point estimate to compare
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vectReps A vector (>=2) of replicate numbers that should be considered as tech reps
  #' @param condName An optional parameter; one-word name for condition
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @param fullOUT Set true if you want full output with all computationally-internal dfs
  #' @return A table of gene names, AIs + CIs, classification into genes demonstrating differential from point estimate AI and those that don't
  #' @examples
  #'

  # Bonferroni correction:
    Q <- 1 - (1-Q)/nrow(inDF)

  if (!fullOUT){
    QCI <- PerformCIAIAnalysis(inDF=inDF, vectReps=vectReps, condName=condName,
                               Q=Q, EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType,
                               fullOUT=F)
  } else {
    OUT <- PerformCIAIAnalysis(inDF=inDF, vectReps=vectReps, condName=condName,
                               Q=Q, EPS=EPS, thr=thr, thrUP=thrUP, thrType=thrType,
                               fullOUT=T)
    QCI <- OUT$AICI
  }

  # Find intersecting intervals > call them FALSE (non-rejected H_0)
  QCI$diffAI <- !(QCI$meanAILow <= pt & QCI$meanAIHigh >= pt)
  QCI$DAE <- (QCI$diffAI & (QCI$meanAI - pt >= minDifference))

  if (!fullOUT){
    return(QCI)
  } else {
    return(list(AICI = QCI,
                FULL_OUT = OUT))
  }
}



