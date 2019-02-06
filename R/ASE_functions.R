# _______________________________________________________________________________________

options(stringsAsFactors = FALSE)
# _______________________________________________________________________________________

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: READ FILES AND CREATE UNIFORM TABLES
# ---------------------------------------------------------------------------------------
#                 TODO: 1. Kallisto functions
# ---------------------------------------------------------------------------------------

GetGatkPipelineTabs <- function(inFiles, nReps, multiple = TRUE, chrom = F){
  #' (GATK pipeline) Concatenate vertically (uniting merge) tables from inFiles.
  #'
  #' @param inFiles A vector of full pathes to files with gene names and allelic counts
  #' @param nReps A vector of numbers, each entry is a number of replicates in the corresponding file
  #' @param multiple Parameter defining if multiple input files are used, default set to TRUE
  #' @param chrom Parameter defining if the resulting table includes chromosome column, default set to FALSE
  #' @return A concatenated table with means/counts, each row corresponds to a gene
  #' @examples
  #'
  # TODO : change naming here
  if (multiple) {
    df <- Reduce(function(x,y){merge(x,y,by="ensembl_gene_id")},
                 lapply(1:length(inFiles), function(x){
                   df0 <- read_delim(inFiles[x], delim="\t", escape_double = FALSE, trim_ws = TRUE)
                   return(df0[,c(1:(2*nReps[x]+1))])
                 }
                 )
    )
    # TODO add chrom option for multiple files
  }
  else {
    df <- read_delim(inFiles, delim="\t", escape_double = FALSE, trim_ws = TRUE)
    df_chrom <- df[,c("ensembl_gene_id","chr")]
    df <- df[,c(1:(2*sum(nReps)+1))]
    df <- as.data.frame(df)
    if (chrom) {
      df <- merge(df, df_chrom, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
    }
  }
  nameColumns <- function(rep_n)  {
    paste0("rep", rep(1:rep_n, each = 2), rep(c("_ref", "_alt"), rep_n))
  }
  if (multiple) {
    names(df) <- c("ensembl_gene_id",
                   paste0("rep", rep(1:sum(nReps), each=2), c("_ref", "_alt")))
    if (chrom) {
      names(df)[ncol(df)] <- "chr"
    }
  }
  else {
    names(df) <- c("ensembl_gene_id", unlist(sapply(nReps, nameColumns)))
    if (chrom) {
      names(df)[ncol(df)] <- "chr"
    }
  }
  return(df)
}

GetGatkPipelineSNPTabs <- function(inFiles, nReps){
  #' (GATK pipeline) Concatenate vertically (uniting merge) tables from inFiles
  #'
  #' @param inFiles A vector of pathes to files
  #' @param nReps A vector of numbers, each -- number of replicates in corresponding file
  #' @return A technical replicates-concatenated table with SNPs, rows corresponds to SNPs
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

CountsToAI <- function(df, reps=NA, meth="meanOfProportions", thr=NA, thrType="each"){
  #' Calculates allelic imbalances from merged counts over given replicates (ai(sum_reps(gene)))
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df)
  #' @param meth An optional parameter for method to use, either sum(m)/sum(p), or sum(m/p) (default = sum(m)/sum(p))
  #' @param thr An optional parameter for a threshold on mean coverage (default = no thr)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return mean(mean(m_1,...,m_6))_SNP / mean(mean(m_1+p_1,...,m_6+p6))_SNP
  #' @examples
  #'
  if (all(is.na(reps))) {
    reps <- 1:((ncol(df)-1)/2)
  }
  cs  <- sort(c(sapply(reps, function(x){c(x*2, x*2+1)}))) # columns numbers
  ddf <- df[, cs] # df with needed columns and without gene names column
  if (is.na(thr)) {
    thr <- 0
    greaterThanThr <- rep(T, nrow(ddf))
  }
  if (thrType == "each"){
    greaterThanThr <- (sapply(1:length(reps),
                             function(x){rowSums(ddf[, c(2*x-1, 2*x)])}) >= thr)
    greaterThanThr[greaterThanThr==FALSE] <- NA
    greaterThanThr <- rowSums(greaterThanThr)/length(reps)
  } else if (thrType == "average"){
    greaterThanThr <- (rowSums(ddf)/length(reps) >= thr)
    greaterThanThr[greaterThanThr==FALSE] <- NA
  }

  if (ncol(ddf) == 2) { # if 1 replicate
    p <- (ddf[, 1]/rowSums(ddf)) * greaterThanThr
  } else {             # if more than 1 replicates
    if (meth == "mergedToProportion") {
      ref <- rowSums(ddf[, seq(1, ncol(ddf), 2)])
      alt <- rowSums(ddf[, seq(2, ncol(ddf), 2)])
      p   <- (ref/(ref + alt)) * greaterThanThr
    } else if (meth == "meanOfProportions") {
      aitab <- do.call(cbind, lapply(1:length(reps), function(i){
        ddf[, i*2-1]/(ddf[, i*2-1]+ddf[, i*2])
      }))
      p <- rowMeans(aitab) * greaterThanThr
    }
  }
  p[is.nan(p)] <- NA
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

MergeSumCounts <- function(df, reps=NA){
  #' Creates a table of sums of maternal and paternal count for given replicates.
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df).
  #' @return Table with sums of mat and pat coverages among given replicates.
  #' @examples
  #'
  if(all(is.na(reps))){
    reps <- 1:((ncol(df)-1)/2)
  }
  if(length(reps)==1){
    return(data.frame(ref_reps = df[, reps*2], alt_reps = df[, reps*2+1]))
  }
  return(data.frame(ref_reps = rowSums(df[, reps*2]),
                    alt_reps = rowSums(df[, reps*2+1])))
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: NAMES REWRITING
# ---------------------------------------------------------------------------------------

NumToDoulbledigitChar <- function(x){
  #' Creates double-digit character from a number or a vector of numbers (up t0 99)
  #'
  #' @param x A number from 0 to 99
  #' @return Double-digit character (or vector) from a number ('x' or '0x')
  #' @examples
  #'
  if(length(x)==1){ x = c(x) }
  doubledigitsChars <- sapply(1:length(x), function(i){
    if(x[i]<=9) { return(paste0('0',x[i])) }
    else { return(as.character(x[i])) }
  })
  return(doubledigitsChars)
}

BuildDesign <- function(experimentNames, techReps){
  #' Creates a design matrix for the experiment
  #'
  #' @param experimentNames Vector with names of the experiments
  #' @param techReps Vector with number of technical replicates in each experiment
  #' @return Dataframe with experiments numbered and numbers of columns
  #' @examples
  #'
  rowsSp <- data.frame(matrix(lapply(1:length(techReps), function(x){(2*sum(techReps[1:x-1])+1):(2*sum(techReps[1:x]))}), nrow=length(techReps), byrow=T),stringsAsFactors=FALSE)
  colnames(rowsSp) <- "replicateCols"
  colExp <- data.frame(matrix(lapply(1:length(techReps), function(x){(sum(techReps[1:x])-techReps[x]+1):(sum(techReps[1:x]))}), nrow=length(techReps), byrow=T),stringsAsFactors=FALSE)
  colnames(colExp) <- "replicateNums"
  designMatrix <- cbind(experimentNames, techReps, rowsSp, colExp)
  colnames(designMatrix) <- c("experimentNames", "techReps", "replicateCols", "replicateNums")
  return(designMatrix)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PAIRED COMPARISONS
# ---------------------------------------------------------------------------------------

CreateDeltaAIForAPairRepsDF <- function(df, reps, thrLOW=0, thrUP=F, thr=NA){
  #' Creates a tab with gene mean coverage and AI deltas for a pair of technical replicates
  #'
  #' @param tab A dataframe of genes/transcripts and parental counts for technical replicates in columns
  #' @param reps A parameter for setting 2 replicates for consideration
  #' @param thrLOW Threshold for min gene coverage (default = 0)
  #' @param thrUP Threshold for max gene coverage (default = F)
  #' @param thr An optinal parameter for threshold on counts
  #' @return A tab with gene mean coverage and AI deltas for a pair of technical replicates
  #' @examples
  #'
  ddf <- data.frame(deltaAI = CountsToAI(df, reps[1], thr) - CountsToAI(df, reps[2], thr),
                    AI1 = CountsToAI(df, reps[1], thr),
                    AI2 = CountsToAI(df, reps[2], thr),
                    MeanCov = MeanCoverage(df, reps))
  ddf <- ddf[ddf$MeanCov >= thrLOW, ]
  if (thrUP) {
    ddf <- ddf[ddf$MeanCov < thrUP, ]
  }
  return(ddf)
}
CreateDeltaAIPairwiseDF <- function(df, thrs=2**c(0:12), thrsSide='both', mlns=F, repnums=F, what="noname", thr=NA){
  #' Creates a table of parwise comparisons for technical replicates in given table
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns
  #' @param thrs An optional vector of thresholds (default = 2**c(0:12))
  #' @param thrsSide An optional parameter of threshold side: 'both' or 'low' (default = 'both')
  #' @param mlns An optionsl binary parameter: if the file contains data of millionr read sampling, FALSE or TRUE (default = F)
  #' @param repnums An optional parameter for a range op replicates for consideration (default = all replicates in df)
  #' @param what A name, is needed if not mlns and no names in list (default = "noname")
  #' @param thr Optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @return A table of parwise comparisons for technical replicates in given table
  #' @examples
  #'
  DFThrNPairsRep <- do.call(rbind, lapply(1:(length(thrs)-1), function(ti){ # [for all threshold bins]
    # 1. set threshold windows and replicates:
    thrs <- 2**c(0:12)
    t <- thrs[ti]
    if (thrsSide=='low') {
      t2 <- Inf
    } else {
      t2 <- thrs[ti+1]
    }
    if(repnums==F) {
      repnums <- 1:(ncol(df)%/%2)
    }
    # 2. take all combinations from given replicates...
    combs <- combn(repnums, 2, function(x){x})
    #  ...and perform paiwise ai comparison:
    ddfThrBin <- do.call(rbind, lapply(1:ncol(combs), function(y){ # [for all replicate combinations]
      ddf          <- CreateDeltaAIForAPairRepsDF(df, combs[,y], thrLOW=t, thrUP=t2, thr=thr)
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

CreateMergedDeltaAIPairwiseDF <- function(df, thrs=2**c(0:12), thrsSide='both', mlns=F, repnums=F, what="noname", thr=NA){
  #' Creates a technical replicates-row-concatenated table for pairwise comparisons of technical replicates in given table
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns
  #' @param thrs An optional vector of thresholds (default = 2**c(0:12))
  #' @param thrsSide An optional parameter of threshold side: 'both' or 'low' (default = 'both')
  #' @param mlns An optionsl binary parameter: if the file contains data of millionr read sampling, FALSE or TRUE (default = F)
  #' @param repnums An optional parameter for a range op replicates for consideration, it can be subreplicates in case of files with millions, as far there 5x columns per each (default = all replicates in df)
  #' @param what A name, is needed if not mlns and no names in list (default = "noname")
  #' @param thr Optional parameter; threshold on the overall number of counts (in all replicates combined) for a gene to be considered
  #' @return A technical replicates-row-concatenated table of parwise comparisons for technical replicates.
  #' @examples
  #'
  if (!mlns) { # if it's not list of millions tabs:
    ddf <- CreateDeltaAIPairwiseDF(df, thrs, thrsSide, mlns, repnums, what, thr)
  } else { # if it's list of tabs (mlns):
    subtabs <- names(df)
    ddf <- do.call(rbind, lapply(subtabs, function(s){
      subdf <- data.frame(df[[s]])
      what  <- s
      return(CreateDeltaAIPairwiseDF(subdf, thrs, thrsSide, mlns, repnums, what, thr))
    }))
  }
  return(ddf)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: BINNING AND QUARTILLING
# ---------------------------------------------------------------------------------------

CreateObservedQuantilesDF <- function(df, P, ep, logbase=T, coverageLimit, group=''){
  #' Creates a table with quantiles and numbers of bins for technical replicates for a given table, binned into log intervals
  #'
  #' @param df A dataframe - output of CreateMergedDeltaAIPairwiseDF()
  #' @param P A vector of %-quartiles
  #' @param logbase The binary parameter if we deal with log base (default = T)
  #' @param ep The log base for binning (0^b, 1^b, ...)
  #' @param coverageLimit Gene coverage limit for consideration
  #' @param group An optional name (default = '')
  #' @return A table with quantiles and numbers of bins for technical replicates for a given table
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
    ddfP <- data.frame(coverageBin = covIntervalsStarts,
                       deltaAI = sapply(1:lenCIS, function(i){
                         dai <- df[df$MeanCov >= covIntervals[i] &
                                     df$MeanCov <  covIntervals[i+1], ]$deltaAI
                         if(p == 0){
                           p = 'sd'
                           sd(dai)
                         } else{
                           quantile(dai, p, na.rm = T)
                         }
                       }),
                       binNObservations = sapply(1:lenCIS, function(i){
                         nrow(df[df$MeanCov >= covIntervals[i] &
                                   df$MeanCov <  covIntervals[i+1], ])
                       }),
                       Q = p)
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
FitLmIntercept <- function(inDf, binNObs, morethan = 10, logoutput = TRUE){
  #' Fits linear model to logarithmic data and outputs intercept for the model with slope=1/2 restriction
  #
  #' @param inDf A dataframe - output of CreateObservedQuantilesDF()
  #' @param N_obs_bin Threshold on number of observations per bin
  #' @param morethan Theshold on gene coverage for lm (default = 10)
  #' @param logoutput Return log intercept? (default = true)
  #' @return lm intercept or log2(lm intercept)
  #' @examples
  #'
  df <- inDf[, c('coverageBin','deltaAI','binNObservations')]
  df <- df[df$coverageBin > morethan &
           df$binNObservations > binNObs, ]
  df[, c('coverageBin','deltaAI')] <- log2(df[, c('coverageBin','deltaAI')])
  loglm <- lm(data = df, deltaAI ~ offset(-0.5*coverageBin), na.action=na.exclude)$coefficients
  if(logoutput){
    return(loglm[1])
  } else {
    return(2**loglm[1])
  }
}

#' FitLmIntercept <- function(inDf, binNObs, morethan=10, logoutput=T){
#'   #' Fits linear model to logarithmic data and counts intersept for model with *1/2 restriction
#'   #'
#'   #' @param inDf A dataframe-output of CreateObservedQuantilesDF()
#'   #' @param N_obs_bin Threshold on number of observations per bin
#'   #' @param morethan Theshold on gene coverage for lm (default = 10)
#'   #' @param logoutput Return log intercept? (default = true)
#'   #' @return lm intercept or log2(lm intercept)
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


# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: AI CI ESTIMATES
# ---------------------------------------------------------------------------------------
CreatePMforAI <- function(dfInt, dfAI, dfCov){
  #' Input: two data frames, one including replicate-gene coverages, one with constants for each technical replicates pair
  #'
  #' @param dfInt A table with column "linInt" of correction constants for each replicates combination (rows), the order of rows should be consistent with columns in dfCov, s.t. pairs are alphabetically ordered
  #' @param dfAI A table with columns for each technical replicate, the rows correspond to genes, the values are AI
  #' @param dfCov A table with columns for each technical replicate, the rows correspond to genes, the values are coverage
  #' @return Plus-minus intervals to determine AI Confidence Intervals for each gene
  #' @examples
  #'
  covSumsCombs <- combn(1:ncol(dfCov), 2, function(x){rowSums(dfCov[, x])})
  covSumsCombs[rowMeans(is.na(dfAI))>0, ] <- NA
  invertCovSumsCombs <- 1 / covSumsCombs
  if(ncol(dfCov) == 2){
    qres = apply(invertCovSumsCombs, 1, function(c){c * dfInt$linInt**2})
  } else if(ncol(dfCov) > 2){
    qres = rowSums(t(apply(invertCovSumsCombs, 1, function(c){c * dfInt$linInt**2})))
  }
  return(0.5/(ncol(dfCov)*(ncol(dfCov)-1))*sqrt(qres))
}
CreateCIforAI <- function(dfInt, dfCounts, thr=NA){
  #' Input: two data frames, one including replicate-gene mat|pat coverages, one with constants for each technical replicates pair
  #'
  #' @param dfInt A table with column "linInt" of correction constants for each replicates combination (rows), the order of rows should be consistent with columns in dfCounts, s.t. pairs are alphabetically ordered
  #' @param dfAI A table with columns for each technical replicate, the rows correspond to genes, the values are AI
  #' @param dfCounts A table with pairs of columns for each technical replicate, the rows correspond to genes, the values are mat|pat coverages
  #' @return 5-column df: ID, AI, Plus-minus intervals, AI Confidence Intervals left and right ends
  #' @examples
  #'
  dfAI  <- sapply(1:(ncol(dfCounts)%/%2), function(i){
    CountsToAI(dfCounts, reps=i, thr=thr)
  })
  dfCov <- sapply(1:(ncol(dfCounts)%/%2), function(i){
    dfCounts[, (2*i)] + dfCounts[, (2*i+1)]
  })

  df <- data.frame(
    ID      = dfCounts[,1],
    meanCov = rowMeans(dfCov),
    meanAI  = CountsToAI(dfCounts, meth="meanOfProportions", thr=thr),
    pm      = CreatePMforAI(dfInt, dfAI, dfCov)
  )
  # if(ncol(dfAI) ==  2){
  #   AImin <- apply(dfAI, 1, mean)
  #   AImax <- AImin
  # } else {
  #   AImin <- apply(dfAI, 1, function(a){rev(sort(a))[2]})
  #   AImax <- apply(dfAI, 1, function(a){sort(a)[2]})
  # }
  AImin <- apply(dfAI, 1, mean)
  AImax <- AImin
  df$meanAILow  <- sapply(AImin-df$pm, function(x){max(0, x)})
  df$meanAIHigh <- sapply(AImax+df$pm, function(x){min(1, x)})
  return(df)
}
# .......................................................................................
# Aftercomments:
# _______________________________________________________________________________________


# THE END
