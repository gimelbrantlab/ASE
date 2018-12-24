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
source("ASE_functions.R")
# _______________________________________________________________________________________

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PERFORM DIFF AI ANALYSIS ON 2 CONDITIONS
# ---------------------------------------------------------------------------------------

PerformDiffAIAnalysisFor2Conditions <- function(inDF, vect1CondReps, vect2CondReps, cond1Name="Condition1", cond2Name="Condition2", Q=0.95, thr=NA, fullOUT=F){
  #' Input: data frame with gene names and counts (reference and alternative) + numbers of replicates to use for each condition
  #'
  #' @param inDF A table with ref & alt counts per gene/SNP for each replicate plus the first column with gene/SNP names
  #' @param vect1CondReps A vector (>=2) of replicate numbers that should be considered as first condition's tech reps
  #' @param vect2CondReps A vector (>=2) of replicate numbers that should be considered as second condition's tech reps
  #' @param cond1Name An optional parameter; one-word name for condition 1
  #' @param cond2Name An optional parameter; one-word name for condition 2
  #' @param Q An optional parameter; %-quantile (for example 0.95, 0.8, etc)
  #' @param thr An optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @return A table of gene names, AIs + CIs for each condition, classification into genes demonstrating differential AI and those that don't
  #' @examples
  #'

  # Take subtables for 1 and 2 conditions:
  dfCondition <- list(cond1 = inDF[, sort(c(1, vect1CondReps*2, vect1CondReps*2+1))],
                      cond2 = inDF[, sort(c(1, vect2CondReps*2, vect2CondReps*2+1))])

  # Create pairvise AI differences for all techreps pairs:
  deltaAIPairwiseDF <- rbind(CreateMergedDeltaAIPairwiseDF(dfCondition$cond1, what=cond1Name, thr=thr),
                             CreateMergedDeltaAIPairwiseDF(dfCondition$cond2, what=cond2Name, thr=thr))
  deltaAIPairwiseDF$group <- paste(deltaAIPairwiseDF$what, deltaAIPairwiseDF$ij)

  # Count quartiles for Mean Coverage bins

  observedQuartilesDF <- do.call(rbind,
    lapply(unique(deltaAIPairwiseDF$group),
           function(gr){
             df  <- deltaAIPairwiseDF[deltaAIPairwiseDF$group == gr, ]
             res <- CreateObservedQuantilesDF(df,
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

  # Calculate AI CIs:

  QCI <- data.frame(ID = inDF[, 1],
                    CreateCIforAI(linIntercepts[[cond1Name]], dfCondition$cond1, thr=thr)[, -1],
                    CreateCIforAI(linIntercepts[[cond2Name]], dfCondition$cond2, thr=thr)[, -1]
  )
  names(QCI) = c("ID",
                 "meanCov1", "meanAI1", "pm1", "meanAI1Low", "meanAI1High",
                 "meanCov2", "meanAI2", "pm2", "meanAI2Low", "meanAI2High")

  # Find intersecting intervals > call them FALSE

  QCI$diffAI <- !(QCI$meanAI1Low < QCI$meanAI2Low & QCI$meanAI1High >= QCI$meanAI2Low |
                  QCI$meanAI1Low >= QCI$meanAI2Low & QCI$meanAI1Low <= QCI$meanAI2High)

  if (!fullOUT){
    return(QCI)
  } else {
    return(list(deltaAIPairwise = deltaAIPairwiseDF,
                observedQuartiles = observedQuartilesDF,
                intercepts = linIntercepts,
                AICI = QCI))
  }

}




# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# EXAMPLE TEST:
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
removeX <- function(DF, legitim_chrgenes){
  return(DF[DF$ensembl_gene_id %in% legitim_chrgenes$gene, ])
}
chrgenes = read.delim('../../../data/Mus_musculus.GRCm38.68.chrgenes.txt', col.names = c('chr', 'gene'))

inTabs = paste0("../../../data/full/",
                c("NEB", "SMARTseq10ng", "SMARTseq100pg"),
                "_processed_gene_extended2.txt")
inTab = "../../../data/5aza/pr_20180714_ISEKI_processed_gene_extended2.txt"

inDF18 = removeX(GetGatkPipelineTabs(inTabs, c(6,6,6)), chrgenes)
inDF5aza = removeX(GetGatkPipelineTabs(inTab, c(13), multiple = F), chrgenes)

RESULT18 = PerformDiffAIAnalysisFor2Conditions(inDF18, vect1CondReps=2:3, vect2CondReps=7:9, Q=0.95)
RESULT5aza = PerformDiffAIAnalysisFor2Conditions(inDF5aza, vect1CondReps=3:4, vect2CondReps=7:9, Q=0.95)

RESULT18; RESULT5aza


