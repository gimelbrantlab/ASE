#
# PERFORM DIFF AI ANALYSIS ON 2 CONDITIONS
# _______________________________________________________________________________________
# TODO
#'
#'* add correction of CI on mean(AI_gene)
#'* add correction on overdispersion
#'* add bins->limits
#'* think about renaming columns to rep1_ref and so on, it's bad
#'* NA on 0/0
#'* NA when counts < thr
# _______________________________________________________________________________________

options(stringsAsFactors = FALSE)
source("yet_another_functions_up_to_date_24112018.R")
# _______________________________________________________________________________________

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PERFORM DIFF AI ANALYSIS ON 2 CONDITIONS
# ---------------------------------------------------------------------------------------

PerformDiffAIAnalysisFor2Conditions <- function(inDF, vect1CondReps, vect2CondReps, cond1Name="Condition1", cond2Name="Condition2", Q=0.95, thr=NA){
  #' Input: data frame with gene names and counts (reference and alternative) + number of replicates for each condition
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vect1CondReps A vector (>=2) of replicate numbers that should be considered as first condition's tech reps
  #' @param vect2CondReps A vector (>=2) of replicate numbers that should be considered as second condition's tech reps
  #' @param cond1Name A one-word name for condition 1
  #' @param cond2Name A one-word name for condition 2
  #' @param Q %-quantile (for example 0.95, 0.8, etc)
  #' @param thr Optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @return A table of gene names, AIs + CIs for each condition, classification into genes demonstrating differential AI and those that don't
  #' @examples
  #'

  # Take subtables for 1 and 2 conditions:
  dfCondition <- list(cond1 = inDF[, sort(c(1, vect1CondReps*2, vect1CondReps*2+1))],
                      cond2 = inDF[, sort(c(1, vect2CondReps*2, vect2CondReps*2+1))])

  # Create pairvise AI differences for all techreps pairs:
  deltaAIPairwiseDF <- rbind(CreateMergedDeltaAIPairwiseDF(dfCondition$cond1, what=cond1Name, thr),
                             CreateMergedDeltaAIPairwiseDF(dfCondition$cond2, what=cond2Name, thr))
  deltaAIPairwiseDF$group <- paste(deltaAIPairwiseDF$what, deltaAIPairwiseDF$ij)

  # Count quartiles for Mean Coverage bins

  observedQuartilesDF <- do.call(rbind,
    lapply(unique(deltaAIPairwiseDF$group),
           function(gr){
             df  <- deltaAIPairwiseDF[deltaAIPairwiseDF$group == gr, ]
             res <- CreateObservedQuartilesDF(df,
                                              P=Q, ep=1.3, logbase=T,
                                              coverageLimit=quantile(deltaAIPairwiseDF$MeanCov, 0.995),
                                              group=gr)
           }
    )
  )
  observedQuartilesDF$condition <- sapply(as.character(observedQuartilesDF$group),
                                          function(x){unlist(strsplit(x, ' '))[1]})
  observedQuartilesDF$ij <- sapply(as.character(observedQuartilesDF$group),
                                   function(x){paste(unlist(strsplit(x, ' '))[2:4], collapse=' ')})

  # Count intercepts:
  linIntercepts <- lapply(unique(observedQuartilesDF$condition),
                          function(cond){
     df <- observedQuartilesDF[observedQuartilesDF$condition == cond, ] # кондишны
     res <- sapply(unique(df$ij), function(x){ # проход по всем парам
       FitLmIntercept(df[df$ij == x, ], binNObs=30, morethan=10, logoutput=F)
     })
     data.frame(condition = cond, ij = unique(df$ij),
                linInt = as.double(res))
   })
  names(linIntercepts) <- unique(observedQuartilesDF$condition)
  # -------------------------------------------------------------------------------------
  # Concerns: Meancoverage is not accurate measure.
  #           Overdispersion itself should be enough, in case of 5aza
  #           it's even worse because the different level of coverage.
  # -------------------------------------------------------------------------------------

  # Count %-CI:
  dfAICondition <- list(cond1 = do.call(cbind, lapply(1:length(vect1CondReps),
                                                          function(i){CountsToAI(dfCondition$cond1, reps=i)})),
                        cond2 = do.call(cbind, lapply(1:length(vect2CondReps),
                                                          function(i){CountsToAI(dfCondition$cond2, reps=i)})))
  dfCovCondition <- list(cond1 = do.call(cbind, lapply(1:length(vect1CondReps),
                                                           function(i){rowSums(dfCondition$cond1[,(2*i):(2*i+1)])})),
                         cond2 = do.call(cbind, lapply(1:length(vect2CondReps),
                                                           function(i){rowSums(dfCondition$cond2[,(2*i):(2*i+1)])})))
  aiCondition <- lapply(dfAICondition, rowMeans)

  CreatePMforAI <- function(dfInt, dfCov){
    covSumsCombs <- combn(1:ncol(dfCov), 2, function(x){rowSums(dfCov[, x])})
    invertCovSumsCombs = 1 / covSumsCombs
    if(ncol(dfCov) == 2){
      qres = apply(invertCovSumsCombs, 1, function(c){c * dfInt$linInt**2})
    } else if(ncol(dfCov) > 2){
      qres = rowSums(t(apply(invertCovSumsCombs, 1, function(c){c * dfInt$linInt**2})))
    }
    return(0.5/(ncol(dfCov)*(ncol(dfCov)-1))*sqrt(qres))
  }

  QCI <- data.frame(ID = inDF[, 1],
                    meanAI1 = aiCondition$cond1,
                    pm1 = CreatePMforAI(linIntercepts[[cond1Name]], dfCovCondition$cond1),
                    meanAI2 = aiCondition$cond2,
                    pm2 = CreatePMforAI(linIntercepts[[cond2Name]], dfCovCondition$cond2))

  QCI$pm1[is.infinite(QCI$pm1)] = 1
  QCI$pm2[is.infinite(QCI$pm2)] = 1

  QCI$meanAI1Low  <- sapply(QCI$meanAI1 - QCI$pm1, function(x){max(0, x)})
  QCI$meanAI1High <- sapply(QCI$meanAI1 + QCI$pm1, function(x){min(1, x)})
  QCI$meanAI2Low  <- sapply(QCI$meanAI2 - QCI$pm2, function(x){max(0, x)})
  QCI$meanAI2High <- sapply(QCI$meanAI2 + QCI$pm2, function(x){min(1, x)})
  QCI[, c("meanAI1Low", "meanAI1High", "meanAI2Low", "meanAI2High")]

  QCI$diffAI <- !(QCI$meanAI1Low < QCI$meanAI2Low & QCI$meanAI1High >= QCI$meanAI2Low |
                  QCI$meanAI1Low >= QCI$meanAI2Low & QCI$meanAI1Low <= QCI$meanAI2High)

  QCI[, c("meanAI1Low", "meanAI1High", "meanAI2Low", "meanAI2High")]


  return(QCI)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# EXAMPLE TEST:
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# inTabs = paste0("../../../data/full/",
#                  c("NEB", "SMARTseq10ng", "SMARTseq100pg"),
#                  "_processed_gene_extended2.txt")
# inDF = GetGatkPipelineTabs(inTabs, c(6,6,6))
# RESULT = PerformDiffAIAnalysisFor2Conditions(inDF, vect1CondReps=2:4, vect2CondReps=7:8, Q=0.95)
# head(RESULT)


