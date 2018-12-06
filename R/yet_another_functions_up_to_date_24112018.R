# 
# YET ANOTHER ATTEMPT TO MAKE NORMAL VERSION OF UP-TO-DATE FUNCTIONS 
# _______________________________________________________________________________________

options(stringsAsFactors = FALSE)
# _______________________________________________________________________________________

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: READ FILES AND CREATE UNIFORM TABLES
# ---------------------------------------------------------------------------------------
#                 TODO: 1. Kallisto functions
# ---------------------------------------------------------------------------------------

GetGatkPipelineTabs <- function(inFiles, nReps){
  #' (GATK pipeline) Concatenate vertically (uniting merge) tables from inFiles.
  #' 
  #' @param inFiles A vector of pathes to files.
  #' @param nReps A vector of numbers, each -- number of replicates in corresponding file.
  #' @return A techreps-concatenated table with means/counts, rows corresponds to genes.
  #' @examples
  #' 
  df <- Reduce(function(x,y){merge(x,y,by="ensembl_gene_id")}, 
               lapply(1:length(inFiles), function(x){
                   df0 <- read_delim(inFiles[x], delim="\t", escape_double = FALSE, trim_ws = TRUE)
                   return(df0[,c(1:(2*nReps[x]+1))])
                 }
               )
               ) 
  names(df) <- c("ensembl_gene_id", 
                 paste0("rep", rep(1:sum(nReps), each=2), c("_ref", "_alt")))
  return(df)
}
GetGatkPipelineSNPTabs <- function(inFiles, nReps){
  #' (GATK pipeline) Concatenate vertically (uniting merge) tables from inFiles.
  #' 
  #' @param inFiles A vector of pathes to files.
  #' @param nReps A vector of numbers, each -- number of replicates in corresponding file.
  #' @return A techreps-concatenated table with SNPs, rows corresponds to SNPs.
  #' @examples
  #' 
  df <- Reduce(function(x,y){merge(x,y,by="ensembl_gene_id")}, 
               lapply(1:length(inFiles), function(x){
                   df <- read_delim(inFiles[x], delim="\t", escape_double = FALSE, trim_ws = TRUE)
                   cs <- which(sapply(names(df), function(x){
                                  (grepl("rep", x, fixed=TRUE) & grepl("aggr", x, fixed=TRUE))
                                }
                              )
                         )
                   snpdf <- apply(df[, cs], 2, function(c){as.numeric(unlist(strsplit(c, ', ')))})
                   snpdf <- data.frame(ensembl_gene_id = paste0("SNP", 1:nrow(snptab_list[[1]])), 
                                       snpdf)
                   return(snpdf)
                 }
               )
               ) 
  names(df) <- c("ensembl_gene_id", 
                 paste0("rep", rep(1:sum(nReps), each=2), c("_ref", "_alt")))
  return(df)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: ALLELIC IMBALANSE AND MEAN COVERAGE
# ---------------------------------------------------------------------------------------

CountsToAI <- function(df, reps=NA){
  #' Calculates allelic imbalance from merged counts over given replicates (ai(sum_reps(gene))).
  #' 
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df).
  #' @return mean(mean(m_1,...,m_6))_SNP / mean(mean(m_1+p_1,...,m_6+p6))_SNP.
  #' @examples
  #' 
  if(all(is.na(reps))){
    reps <- 1:((ncol(df)-1)/2)
  }
  cs  <- sort(c(sapply(reps, function(x){c(x*2, x*2+1)})))
  ddf <- df[, cs] # df with needed columns and without gene names col
  
  if(ncol(ddf) == 2) { # if 1 replicate
    ref <- ddf[, 1]
    alt <- ddf[, 2]
  } else {             # if more than 1 replicates
    ref <- rowSums(ddf[, seq(1, ncol(ddf), 2)])
    alt <- rowSums(ddf[, seq(2, ncol(ddf), 2)])
  }
  p <- (ref/(ref + alt))
  p[is.nan(p)] <- 0
  return(p)
}
MeanCoverage <- function(df, reps=NA){
  #' Calculates mean coverage mat+pat among given replicates.
  #' 
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df).
  #' @return Mean coverage mat+pat among given replicates.
  #' @examples
  #' 
  if(all(is.na(reps))){
    reps <- 1:((ncol(df)-1)/2)
  }
  cs <- sort(c(sapply(reps, function(x){c(x*2-1, x*2)})) + 1)
  return(rowMeans(df[, cs], na.rm=T)*2)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: NAMES REWRITING
# ---------------------------------------------------------------------------------------

NumToDoulbledigitChar <- function(x){
  #' Creates double-digit character from a number or a vector of numbers (up t0 99).
  #' 
  #' @param x A number from 0 to 99.
  #' @return Double-digit character (or vector) from a number ('x' or '0x').
  #' @examples
  #' 
  if(length(x)==1){ x = c(x) }
  doubledigitsChars <- sapply(1:length(x), function(i){
    if(x[i]<=9) { return(paste0('0',x[i])) } 
    else { return(as.character(x[i])) }
  })
  return(doubledigitsChars)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PAIRED COMPARISONS
# ---------------------------------------------------------------------------------------

CreateDeltaAIForAPairRepsDF <- function(df, reps, thr=0, thrUP=F){
  #' Creates a tab with gene mean coverage and AI deltas for a pair of technical replicates.
  #' 
  #' @param tab A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param reps A parameter for setting 2 replicates for consideration.
  #' @param thr Threshold for min gene coverage (default = 0).
  #' @param thrUP Threshold for max gene coverage (default = F).
  #' @return A tab with gene mean coverage and AI deltas for a pair of technical replicates.
  #' @examples 
  #' 
  ddf <- data.frame(deltaAI = CountsToAI(df, reps[1]) - CountsToAI(df, reps[2]), 
                    AI1 = CountsToAI(df, reps[1]),
                    AI2 = CountsToAI(df, reps[2]),
                    MeanCov = MeanCoverage(df, reps))
  ddf <- ddf[ddf$MeanCov >= thr, ]
  if(thrUP) {
    ddf <- ddf[ddf$MeanCov < thrUP, ]
  }
  return(ddf)
}
CreateDeltaAIPairwiseDF <- function(df, thrs=2**c(0:12), thrsSide='both', mlns=F, repcolmns=F, what="noname"){
  #' Creates a table of parwise comparisons for techreps in given table.
  #' 
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param thrs An optional vector of thresholds (default = 2**c(0:12)).
  #' @param thrsSide An optional parameter of threshold side: 'both' or 'low' (default = 'both').
  #' @param mlns An optionsl binary parameter: if the file contains data of millionr read sampling, FALSE or TRUE (default = F).
  #' @param repcolmns An optional parameter for a range op replicates for consideration (default = all replicates in df).
  #' @param what A name, is needed if not mlns and no names in list (default = "noname").
  #' @return A table of parwise comparisons for techreps in given table.
  #' @examples
  #' 
  DFThrNPairsRep <- do.call(rbind, lapply(1:(length(thrs)-1), function(ti){ # [for all threshold bins]
    # 1. set threshold windows and replicates:
    t <- thrs[ti] 
    if(thrsSide=='low') {
      t2 <- Inf
    } else {
      t2 <- thrs[ti+1]
    }
    if(repcolmns==F) {
      repcolmns <- 1:(ncol(df)%/%2)
    }
    # 2. take all combinations from given replicates...
    combs <- combn(repcolmns, 2, function(x){x})
    #  ...and perform paiwise ai comparison:
    ddfThrBin <- do.call(rbind, lapply(1:ncol(combs), function(y){ # [for all replicate combinations]
      ddf          <- CreateDeltaAIForAPairRepsDF(df, combs[,y], t, t2)
      ddf$deltaAI  <- abs(ddf$deltaAI)
      ddf$what     <- what
      ddf$iLocal   <- combs[1,y]
      ddf$jLocal   <- combs[2,y]
      ddf$ijLocal  <- paste(NumToDoulbledigitChar(combs[1,y]), 
                            'vs',
                            NumToDoulbledigitChar(combs[2,y]))
      if(mlns==F) {
        ddf$i      <- combs[1,y]
        ddf$j      <- combs[2,y]
        ddf$ij       <- paste(NumToDoulbledigitChar(combs[1,y]), 
                              'vs',
                              NumToDoulbledigitChar(combs[2,y]))
      } else {
        ddf$i      <- combs[1,y]%/%5 + 1
        ddf$j      <- combs[2,y]%/%5 + 1
        ddf$ij       <- paste(NumToDoulbledigitChar(combs[1,y]%/%5 + 1), 
                              'vs',
                              NumToDoulbledigitChar(combs[2,y]%/%5 + 1))
      }
      ddf$thrLow  <- t
      ddf$thrHigh <- t2
      return(ddf)
    }))
    return(ddfThrBin)
  }))
  return(DFThrNPairsRep)
}
CreateMergedDeltaAIPairwiseDF <- function(df, thrs=2**c(0:12), thrsSide='both', mlns=F, repcolmns=F, what="noname"){
  #' Creates a techreps-row-concatenated table of parwise comparisons for techreps in given table, via CreateDeltaAIPairwiseDF().
  #' 
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param thrs An optional vector of thresholds (default = 2**c(0:12)).
  #' @param thrsSide An optional parameter of threshold side: 'both' or 'low' (default = 'both').
  #' @param mlns An optionsl binary parameter: if the file contains data of millionr read sampling, FALSE or TRUE (default = F).
  #' @param repcolmns An optional parameter for a range op replicates for consideration, it can be subreplicates in case of files with millions, as far there 5x columns per each (default = all replicates in df).
  #' @param what A name, is needed if not mlns and no names in list (default = "noname").
  #' @return A techreps-row-concatenated table of parwise comparisons for techreps.
  #' @examples
  #' 
  if (!mlns){ # if it's not list of millions tabs:
    ddf <- CreateDeltaAIPairwiseDF(df, thrs, thrsSide, mlns, repcolmns, what)
  } else { # if it's list of tabs (mlns):
    subtabs <- names(df)
    ddf <- do.call(rbind, lapply(subtabs, function(s){
      subdf <- data.frame(df[[s]])
      what  <- s
      return(CreateDeltaAIPairwiseDF(subdf, thrs, thrsSide, mlns, repcolmns, what))
    }))
  }
  return(ddf)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: BINNING AND QUARTILLING
# ---------------------------------------------------------------------------------------

CreateObservedQuartilesDF <- function(df, P, ep, logbase=T, coverageLimit, group=''){
  #' Creates a table of parwise AI differences quartilles info for techreps in given table, binned to log intervals.
  #' 
  #' @param df A dataframe-output of CreateMergedDeltaAIPairwiseDF().
  #' @param P A vector of %-quartiles.
  #' @param logbase The binary parameter if we deal with log base (default = T).
  #' @param ep The log base for binning (0^b, 1^b, ...).
  #' @param coverageLimit Gene coverage limit for consideration.
  #' @param group An ptional distinguishible name (default = '').
  #' @return A table of parwise AI differences info for techreps in given table, binned to log intervals.
  #' @examples
  #' 
  ## P=0 is SD
  if(logbase){
    covIntervalsStarts <- unique(floor(ep**(0:log(coverageLimit, base=ep))))
  } else {
    covIntervalsStarts <- seq(0, coverageLimit-ep, ep)
  }
  lenCIS       <- length(covIntervalsStarts)
  covIntervals <- c(covIntervalsStarts, coverageLimit)
  
  
  ddf <- do.call(rbind, lapply(P, function(p){ # [for all quartile (%)]:
    # [for all coverage bins]:
    ddfP <- data.frame(coverageBin     = covIntervalsStarts,
                       deltaAI         = sapply(1:lenCIS, function(i){
                         dai <- df[df$MeanCov >= covIntervals[i] &
                                   df$MeanCov <  covIntervals[i+1], ]$deltaAI
                         if(p == 0){
                           p = 'sd'
                           sd(dai)
                         } else{
                           quantile(dai, p) 
                         }
                       }),
                       binNObservations = sapply(1:lenCIS, function(i){
                         nrow(df[df$MeanCov >= covIntervals[i] &
                                 df$MeanCov <  covIntervals[i+1], ])
                       }),
                       Q                = p)
    return(ddfP)
  })
  )
  ddf$group <- as.factor(paste(group, ddf$Q))
  ddf$Q     <- as.factor(ddf$Q)
  return(ddf)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: FITTING LM
# ---------------------------------------------------------------------------------------

# Was: fit_lm_intercept_how_we_want_morethan(inDf, N_obs_bin, morethan=10)
FitLmIntercept <- function(inDf, binNObs, morethan=10, logoutput=T){
  #' Fits linear model to logarithmic data and counts intersept for model with *1/2 restriction.
  #'
  #' @param inDf A dataframe-output of CreateObservedQuartilesDF().
  #' @param N_obs_bin Threshold on number of observations per bin.
  #' @param morethan Theshold on gene coverage for lm (default = 10).
  #' @param logoutput Return log intercept? (default = true)
  #' @return lm intercept or log2(lm intercept).
  #' @examples
  #'
  df <- inDf[, c('coverageBin','deltaAI','binNObservations')]
  df <- df[df$coverageBin > morethan &
           df$binNObservations > binNObs, ]
  df[, c('coverageBin','deltaAI')] <- log2(df[, c('coverageBin','deltaAI')])
  loglm <- lm(data = df, deltaAI ~ offset(-0.5*coverageBin))$coefficients
  if(logoutput){
    return(loglm[1])
  } else {
    return(2**loglm[1])
  }
}

#' FitLmIntercept <- function(inDf, binNObs, morethan=10, logoutput=T){
#'   #' Fits linear model to logarithmic data and counts intersept for model with *1/2 restriction. 
#'   #' 
#'   #' @param inDf A dataframe-output of CreateObservedQuartilesDF().
#'   #' @param N_obs_bin Threshold on number of observations per bin.
#'   #' @param morethan Theshold on gene coverage for lm (default = 10).
#'   #' @param logoutput Return log intercept? (default = true)
#'   #' @return lm intercept or log2(lm intercept).
#'   #' @examples
#'   #' 
#'   df <- inDf[, c('coverageBin','deltaAI','binNObservations')]
#'   df <- df[df$coverageBin > morethan & 
#'              df$binNObservations > binNObs, ]
#'   df[, c('coverageBin','deltaAI')] <- log2(df[, c('coverageBin','deltaAI')])
#'   # if(sum(is.na(df[, c('coverageBin','deltaAI')]))>0){
#'   #   return(NA)
#'   # }
#'   loglm <- lm(data = df, deltaAI ~ offset(-0.5*coverageBin))$coefficients
#'   if(logoutput){ 
#'     return(loglm[1]) 
#'   } else {
#'     return(2**loglm[1])
#'   }
#' }

# .......................................................................................
# Aftercomments:
# _______________________________________________________________________________________


# THE END