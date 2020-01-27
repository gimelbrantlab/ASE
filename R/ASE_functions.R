#
# BASIC FUNCTIONS
# _______________________________________________________________________________________

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: READ FILES AND CREATE UNIFORM TABLES
# ---------------------------------------------------------------------------------------

GetGatkPipelineTabs <- function(inFiles, nReps, contigs = vector()){
  #' (GATK pipeline) Concatenate vertically (uniting merge) tables from inFiles.
  #'
  #' @param inFiles A vector of full pathes to files with gene names and allelic counts
  #' @param nReps A vector of numbers, each entry is a number of replicates in the corresponding file
  #' @param contigs Parameter defining if the resulting table should be filtered by contig column preserving only row corresponding to a given vector, default set to empty vector() and no filtering applied
  #' @return A concatenated table with means/counts, each row corresponds to a gene
  #' @examples
  #'
  #'
  options(stringsAsFactors = FALSE)
  nameColumns <- function(nReps)  {
    if(length(nReps) != 1){
      lapply(1:length(nReps), function(i){
        paste0("rep", i, ".", rep(1:nReps[i], each = 2), c("_ref", "_alt"))
      })
    } else {
      list(paste0("rep", rep(1:sum(nReps), each = 2), c("_ref", "_alt")))
    }
  }
  if(length(contigs) != 0){
    # i.e. filtering by contig is needed
    cs_merge <- c("ID", "contig")
    # contig column identification:
    c_contig <- which(names(read.delim(inFiles[1], sep="\t")) == "contig")
    if(length(inFiles) == 1){
      cs <- list(c(1, 2:(2*sum(nReps)+1), c_contig))
      cs_names <- list(c("ID", unlist(nameColumns(nReps)), "contig"))
    } else {
      cs <- lapply(1:length(nReps), function(i){c(1, 2:(2*nReps[i]+1), c_contig)})
      cs_names <- lapply(nameColumns(nReps), function(x){c("ID", x, "contig")})
    }
  } else {
    # i.e. NO filtering by contig
    cs_merge <- c("ID")
    if(length(inFiles) == 1){
      cs <- list(c(1, 2:(2*sum(nReps)+1)))
      cs_names <- list(c("ID", unlist(nameColumns(nReps))))
    } else {
      cs <- lapply(1:length(nReps), function(i){c(1, 2:(2*nReps[i]+1))})
      cs_names <- lapply(nameColumns(nReps), function(x){c("ID", x)})
    }
  }

  df <- Reduce(function(x,y){merge(x, y, by=cs_merge)},
               lapply(1:length(inFiles), function(i){
                 df0 <- read.delim(inFiles[i], sep="\t")
                 df0 <- df0[, cs[[i]]]
                 names(df0) <- cs_names[[i]]
                 return(df0)
               })
  )
  if(length(contigs) != 0){
    df <- df[df$contig %in% contigs, ]
    df <- df[, -which(names(df)=="contig")]
    df
  }

  return(df)
}


# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: ALLELIC IMBALANSE AND MEAN COVERAGE
#                            (COMPUTING, FILTERING, REORGANIZING)
# ---------------------------------------------------------------------------------------

ThresholdingCounts <- function(df, reps=NA, thr=NA, thrUP=NA, thrType="each"){
  #' Takes table with gene names and counts and returns table, where all genes that under threshold have NA coverage. Can be restricted to particular reps.
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df)
  #' @param thr An optional parameter for a threshold on gene coverage (default = NA)
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return Table with masked with NA undercovered genes
  #' @examples
  #'
  options(stringsAsFactors = FALSE)
  # Taking replicates:
  if(anyNA(reps)){
    reps <- 1:(ncol(df)%/%2)
    ddf  <- df
  } else {
    cs  <- sort(c(sapply(reps, function(x){c(x*2, x*2+1)}))) # columns numbers
    ddf <- df[, c(1, cs)]
  }

  # Thresholding:
  if(!anyNA(thr)){
    if(thrType == "each"){
      # Masking with NA gene coverage info replicate-specific if it is < thr
      greaterThanThr <- (sapply(1:length(reps), function(x){
          rep_thr_info <- rowSums(ddf[, c(2*x, 2*x+1)])
          cbind(rep_thr_info, rep_thr_info)
        }) >= thr)
      greaterThanThr[greaterThanThr==FALSE] <- NA
      ddf[, -1] <- ddf[, -1] * greaterThanThr
    } else if(thrType == "average"){
      # Masking with NA entire gene coverage info for all replicates if average coverage < thr:
      greaterThanThr <- (rowMeans(ddf[, -1])*2 >= thr)
      greaterThanThr[greaterThanThr==FALSE] <- NA
      ddf[, -1] <- ddf[, -1] * greaterThanThr
    }
  }

  #Same for UPPER Threshold:
  if(!anyNA(thrUP)){
    if(thrType == "each"){
      # Masking with NA gene coverage info replicate-specific if it is > thrUP
      lesserThanThrUP <- (sapply(1:length(reps), function(x){
        rep_thr_info <- rowSums(ddf[, c(2*x, 2*x+1)])
        cbind(rep_thr_info, rep_thr_info)
      }) <= thrUP)
      lesserThanThrUP[lesserThanThrUP==FALSE] <- NA
      ddf[, -1] <- ddf[, -1] * lesserThanThrUP
    } else if(thrType == "average"){
      # Masking with NA entire gene coverage info for all replicates if average coverage > thr:
      lesserThanThrUP <- (rowMeans(ddf[, -1])*2 <= thrUP)
      lesserThanThrUP[lesserThanThrUP==FALSE] <- NA
      ddf[, -1] <- ddf[, -1] * lesserThanThrUP
    }
  }

  return(ddf)
}

CountsToAI <- function(df, reps=NA, meth="meanOfProportions", thr=NA, thrUP=NA, thrType="each"){
  #' Calculates allelic imbalances from merged counts over given replicates (ai(sum_reps(gene)))
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df)
  #' @param meth An optional parameter for method to use, either sum(m)/sum(p), or sum(m/p) (default = sum(m)/sum(p))
  #' @param thr An optional parameter for a threshold on gene coverage (default = NA)
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return Df with names and mean(mean(m_1,...,m_6))_SNP / mean(mean(m_1+p_1,...,m_6+p6))_SNP
  #' @examples
  #'
  options(stringsAsFactors = FALSE)

  ddf <- ThresholdingCounts(df, reps, thr, thrUP, thrType)

  if(ncol(ddf) == 3){ # if 1 replicate
    p <- (ddf[, 2]/rowSums(ddf[, -1]))
  } else {             # if more than 1 replicates
    if (meth == "mergedToProportion") {
      ref <- rowSums(ddf[, seq(2, ncol(ddf), 2)])
      all <- rowSums(ddf[, 2:ncol(ddf)])
      p   <- ref/all
    } else if(meth == "meanOfProportions"){
      aitab <- sapply(1:(ncol(ddf)%/%2), function(i){
        ddf[, i*2]/(ddf[, i*2]+ddf[, i*2+1])
      })
      p <- rowMeans(aitab)
    }
  }
  p[is.nan(p)] <- NA

  res_df <- data.frame(df[, 1], AI = p)
  names(res_df)[1] <- names(df)[1]

  return(res_df)
}

MeanCoverage <- function(df, reps=NA, thr=NA, thrUP=NA, thrType="each"){
  #' Calculates mean coverage mat+pat among given replicates.
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df).
  #' @param thr An optional parameter for a threshold on gene coverage (default = NA)
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return Df with names and Mean coverage mat+pat among given replicates.
  #' @examples
  #'
  options(stringsAsFactors = FALSE)

  ddf <- ThresholdingCounts(df, reps, thr, thrUP, thrType)

  res_df <- data.frame(df[, 1], meanCOV = rowMeans(ddf[, -1])*2)
  names(res_df)[1] <- names(df)[1]

  return(res_df)
}

MergeSumCounts <- function(df, reps=NA, thr=NA, thrUP=NA, thrType="each"){
  #' Creates a table of sums of maternal and paternal count for given replicates.
  #'
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df).
  #' @param thr An optional parameter for a threshold on gene coverage (default = NA)
  #' @param thrUP An optional parameter for a threshold for max gene coverage (default = NA)
  #' @param thrType An optional parameter for threshold type (default = "each", also can be "average" coverage on replicates)
  #' @return Table with names and sums of mat and pat coverages among given replicates.
  #' @examples
  #'
  options(stringsAsFactors = FALSE)

  ddf <- ThresholdingCounts(df, reps, thr, thrUP, thrType)

  if(ncol(ddf) == 3){
    res_df <- data.frame(ddf[, 1], ref_reps = ddf[, 2], alt_reps = ddf[, 3])
  } else {
    res_df <- data.frame(ddf[, 1],
                         ref_reps = rowSums(ddf[, seq(2, ncol(ddf), 2)]),
                         alt_reps = rowSums(ddf[, seq(3, ncol(ddf), 2)]))
  }
  names(res_df)[1] <- names(ddf)[1]

  return(res_df)
}


# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: DESIGN MATRIX
# ---------------------------------------------------------------------------------------

BuildDesign <- function(experimentNames, techReps, corrConst=NA){
  #' Creates a design matrix for the experiment
  #'
  #' @param experimentNames Vector with names of the experiments
  #' @param techReps Vector with number of technical replicates in each experiment
  #' @param corrConst Optional, Vector with correction constants for each experiment
  #' @return Dataframe with experiments numbered and numbers of columns
  #' @examples
  #'
  rowsSp <- data.frame(matrix(lapply(1:length(techReps), function(x){
                                                           (2*sum(techReps[1:x-1])+1):(2*sum(techReps[1:x]))
                                                         }),
                              nrow=length(techReps), byrow=T), stringsAsFactors=FALSE)
  colnames(rowsSp) <- "replicateCols"
  colExp <- data.frame(matrix(lapply(1:length(techReps), function(x){
                                                           (sum(techReps[1:x])-techReps[x]+1):(sum(techReps[1:x]))
                                                         }),
                              nrow=length(techReps), byrow=T), stringsAsFactors=FALSE)
  colnames(colExp) <- "replicateNums"
  designMatrix <- cbind(experimentNames, techReps, rowsSp, colExp)
  colnames(designMatrix) <- c("experimentNames", "techReps", "replicateCols", "replicateNums")
  if (!sum(is.na(corrConst))) {
    designMatrix <- cbind(designMatrix, corrConst)
  }
  return(designMatrix)
}

NameColumns <- function(exp_n, rep_n)  {
  #' Helper function to quickly rename columns in geneCountTab dataframe
  #'
  #' @param exp_n Experiment number
  #' @param rep_n Number of replicates for the experiment
  #' @return Vector with names
  #' @examples colnames(geneCountTab)[2:13] <- c(nameColumns(1,2), nameColumns(2,2), nameColumns(3,2))
  #'

  paste0("exp", rep(exp_n, 2*rep_n), "_rep", rep(1:rep_n, each = 2), "_", rep(c("ref", "alt"), rep_n))
}


# _______________________________________________________________________________________
