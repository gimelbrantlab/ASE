# 
# YET ANOTHER ATTEMPT TO MAKE NORMAL VERSION OF UP-TO-DATE FUNCTIONS 
# _______________________________________________________________________________________

# =======================================================================================
#                 0. Libraries
# =======================================================================================
library("tidyverse")
library("gridExtra")
library("optparse")
# _______________________________________________________________________________________


# =======================================================================================
#                 1. Overall parameters
# =======================================================================================
options(stringsAsFactors = FALSE)
# _______________________________________________________________________________________


# =======================================================================================
#                 2. Functions itself
# =======================================================================================
# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: READ FILES AND CREATE UNIFORM TABLES
# ---------------------------------------------------------------------------------------

readTabMeans = function(inFile, nReps){
  #' (GATK pipeline) Reads file _extended2.txt and takes means/counts columns only.
  #' 
  #' @param inFiles A vector of pathes to files.
  #' @param nReps A vector of numbers, each -- number of replicates in corresponding file.
  #' @return A table with means/counts columns only, rows corresponds to genes.
  #' @examples
  #' 
  read_delim(inFile, delim="\t", escape_double = FALSE, trim_ws = TRUE)[,c(1:(2*nReps+1))]
}
getTabs = function(inFiles, nReps){
  #' (GATK pipeline) Concatenate vertically (uniting merge) tables from inFiles, processed via readTabMeans().
  #' 
  #' @param inFiles A vector of pathes to files.
  #' @param nReps A vector of numbers, each -- number of replicates in corresponding file.
  #' @return A techreps-concatenated table with means/counts, rows corresponds to genes.
  #' @examples
  #' 
  tab = Reduce(function(x,y){merge(x,y,by="ensembl_gene_id")}, 
               lapply(1:length(inFiles), function(x){readTabMeans(inFiles[x], nReps[x])})) 
  names(tab) = c("ensembl_gene_id", paste0("rep", rep(1:sum(nReps), each=2), c("_ref", "_alt")))
  return(tab)
}

readTabKallisto = function(inFile){
  #' (Kallisto pipeline) Reads kallisto output file and perform it in tr-mat-pat tab.
  #' 
  #' @param inFiles A vector of pathes to files.
  #' @return A table with parental counts, rows corresponds to transcripts.
  #' @examples
  #' 
  tab0 = read_delim(inFile, delim="\t", escape_double = FALSE, trim_ws = TRUE, col_types = "ciddd")[, c(1,4)]
  tab0$parent = unlist(lapply(strsplit(tab0$target_id, "_"), function(x){x[2]}))
  tab0$transcript = unlist(lapply(strsplit(tab0$target_id, "_"), function(x){x[1]}))
  tab = merge(tab0[tab0$parent=="129S1",c(4,2)], tab0[tab0$parent=="CAST",c(4,2)], by="transcript")
  names(tab) = c("transcript", "ref", "alt")
  return(tab)
}
createKallistoTrMetaTab = function(inFiles){
  #' (Kallisto pipeline) Reads kallisto output files and creates overall tr-mat1-pat1-mat2-... tab, each processed via readTabKallisto().
  #' 
  #' @param inFiles A vector of pathes to files.
  #' @return A techreps-concatenated table with parental counts, rows corresponds to transcripts.
  #' @examples
  #' 
  tab = Reduce(function(x,y){merge(x,y,by="transcript")}, 
               lapply(inFiles, function(x){readTabKallisto(x)})) 
  names(tab) = c("ensembl_transcript_id", paste0("rep", rep(1:length(inFiles), each=2), c("_ref", "_alt")))
  return(tab)
}

#     SNP Tab:
readSNPtab = function(inFile){
  #' (GATK pipeline) Reads file _extended2.txt and takes SNP data. Woks only for our pipeline!!!
  #' 
  #' @param inFiles A vector of pathes to files.
  #' @return A table with SNP columns only, rows corresponds to SNPs.
  #' @examples
  #' 
  df = read_delim(inFile, delim="\t", escape_double = FALSE, trim_ws = TRUE)
  cs = which(sapply(names(df), function(x){(grepl("rep", x, fixed=TRUE) & grepl("aggr", x, fixed=TRUE))}))
  snpdf = apply(df[, cs], 2, function(c){as.numeric(unlist(strsplit(c, ', ')))})
  snpdf = data.frame(ensembl_gene_id = paste0("SNP", 1:nrow(snptab_list[[1]])), snpdf)
  return(snpdf)
}
getSNPTabs = function(inFiles, nReps){
  #' (GATK pipeline) Concatenate vertically (uniting merge) tables from inFiles, processed via readSNPtab().
  #' 
  #' @param inFiles A vector of pathes to files.
  #' @param nReps A vector of numbers, each -- number of replicates in corresponding file.
  #' @return A techreps-concatenated table with SNPs, rows corresponds to SNPs.
  #' @examples
  #' 
  tab = Reduce(function(x,y){merge(x,y,by="ensembl_gene_id")}, 
               lapply(1:length(inFiles), function(x){readSNPtab(inFiles[x])})) 
  names(tab) = c("ensembl_gene_id", paste0("rep", rep(1:sum(nReps), each=2), c("_ref", "_alt")))
  return(tab)
}


# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: ALLELIC IMBALANSE AND MEAN COVERAGE
# ---------------------------------------------------------------------------------------

genemeansToAI = function(df, reps=NA){
  #' Calculates allelic imbalance from merged counts over given replicates (ai(sum_reps(gene))).
  #' 
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df).
  #' @return mean(mean(m_1,...,m_6))_SNP / mean(mean(m_1+p_1,...,m_6+p6))_SNP.
  #' @examples
  #' 
  if(all(is.na(reps))){
    reps = 1:((ncol(df)-1)/2)
  }
  c = sort(c(sapply(reps, function(x){c(x*2,x*2+1)})))
  ddf = df[,c]
  
  if(ncol(ddf) == 2) {
    ref = ddf[, 1]
    alt = ddf[, 2]
  } else {
    ref = rowSums(ddf[,seq(1, ncol(ddf), 2)], na.rm=T)
    alt = rowSums(ddf[,seq(2, ncol(ddf), 2)], na.rm=T)
  }
  p = (ref/(ref + alt))
  p[is.nan(p)] = 0
  return(p)
}
meanCoverage = function(df, reps=NA){
  #' Calculates mean coverage mat+pat among given replicates.
  #' 
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param reps An optional parameter for a range op replicates for consideration (default = all replicates in df).
  #' @return Mean coverage mat+pat among given replicates.
  #' @examples
  #' 
  if(is.na(reps)){
    reps = 1:((ncol(df)-1)/2)
  }
  c = sort(c(sapply(reps, function(x){c(x*2-1,x*2)})) + 1)
  return(rowMeans(df[,c], na.rm=T)*2)
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: NAMES REWRITING
# ---------------------------------------------------------------------------------------

num_to_doulbledigit_char = function(x){
  #' Creates double-digit character from a number (up t0 99).
  #' 
  #' @param x A number from 0 to 99.
  #' @return Double-digit character from a number ('x' or '0x').
  #' @examples
  #' 
  if(x<=9) { return(paste0('0',x)) } 
  else { return(as.character(x)) }
}
num_to_doulbledigit_vectchar = function(x){
  #' Creates a number from double-digit character.
  #' 
  #' @param x A character 'xy' for x,y from 0 to 9.
  #' @return Number from double-digit character.
  #' @examples
  #' 
  return(sapply(x, function(i){num_to_doulbledigit_char(i)}))
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: PAIRED COMPARISONS
# ---------------------------------------------------------------------------------------

meancovAI2ForPairDF = function(tab, reps, thr=0, thrUP=F){
  #' Creates a tab with gene mean coverage and AI deltas for a pair of technical replicates.
  #' 
  #' @param tab A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param reps A parameter for setting 2 replicates for consideration.
  #' @param thr Threshold for min gene coverage.
  #' @param thrUP Threshold for max gene coverage.
  #' @return A tab with gene mean coverage and AI deltas for a pair of technical replicates.
  #' @examples 
  #' 
  DF = data.frame(genemeansToAI(tab, reps[1]) - genemeansToAI(tab, reps[2]), 
                  genemeansToAI(tab, reps[1]),
                  genemeansToAI(tab, reps[2]),
                  meanCoverage(tab, reps))
  names(DF) = c("diff_AI", "AI1", "AI2", "mean_cov")
  DF = DF[DF$mean_cov >= thr, ]
  if(thrUP) {DF = DF[DF$mean_cov < thrUP, ]}
  return(DF)
}

get_one_df_AIpairdist_distr = function(df, thrs=2**c(0:12), thrs_side='both', mlns=F, repcolmns=F, what="noname"){
  #' Creates a table of parwise comparisons for techreps in given table.
  #' 
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param thrs An optional vector of thresholds (default = 2**c(0:12)).
  #' @param thrs_side An optional parameter of threshold side: 'both' or 'low' (default = 'both').
  #' @param mlns An optionsl binary parameter: if the file contains data of millionr read sampling, FALSE or TRUE (default = F).
  #' @param repcolmns An optional parameter for a range op replicates for consideration (default = all replicates in df).
  #' @param what A name, is needed if not mlns and no names in list (default = "noname").
  #' @return A table of parwise comparisons for techreps in given table.
  #' @examples
  #' 
  do.call(rbind, lapply(1:(length(thrs)-1), function(ti){
    # set threshold windows:
    t = thrs[ti] 
    if(thrs_side=='low') {
      t2 = Inf
    } else {
      t2 = thrs[ti+1]
    }
    if(repcolmns==F) {
      repcolmns = 1:(ncol(df)%/%2)
    }
    # all combinations from given replicates:
    combs = combn(repcolmns, 2, function(x){x})
    do.call(rbind, lapply(1:ncol(combs), function(y){
      df = meancovAI2ForPairDF(df, combs[,y], t, t2)
      df$diff_AI = abs(df$diff_AI)
      df$what = what
      df$i_local = combs[1,y]
      df$j_local = combs[2,y]
      df$ij_local = paste(num_to_doulbledigit_char(combs[1,y]), 'vs' ,
                          num_to_doulbledigit_char(combs[2,y]))
      if(mlns==F) {
        df$ij = paste(num_to_doulbledigit_char(combs[1,y]), 'vs' ,
                      num_to_doulbledigit_char(combs[2,y]))
        df$i = combs[1,y]
        df$j = combs[2,y]
      } else {
        df$i = combs[1,y]%/%5 + 1
        df$j = combs[2,y]%/%5 + 1
        df$ij = paste(num_to_doulbledigit_char(combs[1,y]%/%5 + 1), 'vs' ,
                      num_to_doulbledigit_char(combs[2,y]%/%5 + 1))
      }
      df$thr_low = t
      df$thr_high = t2
      return(df)
    }))
  }))
}

do_df_AIpairdist_distr = function(df, thrs=2**c(0:12), thrs_side='both', mlns=F, repcolmns=F, what="noname"){
  #' Creates a techreps-row-concatenated table of parwise comparisons for techreps in given table, via get_one_df_AIpairdist_distr().
  #' 
  #' @param df A dataframe of genes/transcripts and parental counts for technical replicates in columns.
  #' @param thrs An optional vector of thresholds (default = 2**c(0:12)).
  #' @param thrs_side An optional parameter of threshold side: 'both' or 'low' (default = 'both').
  #' @param mlns An optionsl binary parameter: if the file contains data of millionr read sampling, FALSE or TRUE (default = F).
  #' @param repcolmns An optional parameter for a range op replicates for consideration (default = all replicates in df).
  #' @param what A name, is needed if not mlns and no names in list (default = "noname").
  #' @return A techreps-row-concatenated table of parwise comparisons for techreps.
  #' @examples
  #' 
  # if it's not list of millions tabs:
  if (!mlns){
    get_one_df_AIpairdist_distr(df, thrs, thrs_side, mlns, repcolmns, what)
  }
  # if it's list of tabs:
  else {
    subtabs = names(df)
    do.call(rbind, lapply(subtabs, function(s){
      subdf = data.frame(df[[s]])
      what = s
      get_one_df_AIpairdist_distr(subdf, thrs, mlns, thrs_side, repcolmns, what)
    }))
  }
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: BINNING AND QUARTILLING
# ---------------------------------------------------------------------------------------

df_obs_quartiles_info = function(df, P, ep, n){
  #' Creates a table of parwise AI differences quartilles info for techreps in given table, binned to linear intervals.
  #' 
  #' @param df A dataframe-output of do_df_AIpairdist_distr().
  #' @param P A vector of %-quartiles.
  #' @param ep A parameter for gene coverage window for binning.
  #' @param n Gene coverage limit for consideration.
  #' @return A table of parwise AI differences info for techreps in given table, binned to linear intervals.
  #' @examples
  #' 
  df = do.call(rbind,
               lapply(1:length(P), function(i){
                 p = P[i]
                 sdc = qnorm((1-p)/2 + p)
                 data.frame(MeanGeneCoverage = seq(0, n-ep, ep) + ep/2, 
                            deltaAI = sapply(seq(0, n-ep, ep), function(x){
                              dai = df[df$mean_cov >= x & df$mean_cov < x+ep, ]$diff_AI
                              quantile(dai, p*nonout)
                            }), 
                            N_observations = sapply(seq(0, n-ep, ep), function(x){
                              nrow(df[df$mean_cov >= x & df$mean_cov < x+ep, ])
                            }),
                            Q = p, Type = 'observed')
               })
  )
  df$group = as.factor(paste(df$group, df$Type, df$Q))
  df$Q = as.factor(df$Q)
  df$Type = as.factor(df$Type)
  df
}
df_obs_logbase_quartiles_info = function(df, P, logbase, n, group=''){
  #' Creates a table of parwise AI differences quartilles info for techreps in given table, binned to log intervals.
  #' 
  #' @param df A dataframe-output of do_df_AIpairdist_distr().
  #' @param P A vector of %-quartiles.
  #' @param logbase The log base for binning (0^b, 1^b, ...).
  #' @param n Gene coverage limit for consideration.
  #' @param group An ptional distinguishible name (default = '').
  #' @return A table of parwise AI differences info for techreps in given table, binned to log intervals.
  #' @examples
  #' 
  ## P=0 is SD
  df = do.call(rbind,
               lapply(1:length(P), function(i){
                 p = P[i]
                 meancov_intervals_start = unique(floor(logbase**(0:log(n, base=logbase))))
                 lmi = length(meancov_intervals_start)
                 meancov_intervals = c(meancov_intervals_start, n)
                 data.frame(MeanGeneCoverage = meancov_intervals_start,
                            deltaAI = sapply(1:lmi, function(x){
                              dai = df[df$mean_cov >= meancov_intervals[x] &
                                         df$mean_cov < meancov_intervals[x+1], ]$diff_AI
                              if(p == 0){
                                p = 'sd'
                                sd(dai)
                              }
                              else{
                                quantile(dai, p) #p*nonout
                              }
                            }),
                            N_observations = sapply(1:lmi, function(x){
                              nrow(df[df$mean_cov >= meancov_intervals[x] &
                                        df$mean_cov < meancov_intervals[x+1], ])
                            }),
                            Q = p, Type = 'observed')
               })
  )
  df$group = as.factor(paste(group, df$Type, df$Q))
  df$Q = as.factor(df$Q)
  df$Type = as.factor(df$Type)
  df
}

# ---------------------------------------------------------------------------------------
#                 FUNCTIONS: FITTING LM
# ---------------------------------------------------------------------------------------

fit_lm_intercept_how_we_want_morethan = function(in_df, N_obs_bin, morethan=10){
  #' Fits linear model to logarithmic data and counts intersept for model with *1/2 restriction. 
  #' 
  #' @param in_df A dataframe-output of df_obs_logbase_quartiles_info() or df_obs_quartiles_info().
  #' @param N_obs_bin Threshold on number of observations per bin.
  #' @param morethan Theshold on gene coverage for lm (default = 10).
  #' @return lm intercept.
  #' @examples
  #' 
  df = in_df[in_df$MeanGeneCoverage > morethan & in_df$N_observations > N_obs_bin, 
             c('MeanGeneCoverage','deltaAI','N_observations')]
  df[, c('MeanGeneCoverage','deltaAI')] = log2(df[, c('MeanGeneCoverage','deltaAI')])
  2**lm(data = df, deltaAI ~ offset(-0.5*MeanGeneCoverage))$coefficients[1]
}
fit_lmlog_intercept_how_we_want_morethan = function(in_df, N_obs_bin, morethan=10){
  #' Fits linear model to logarithmic data and counts intersept for model with *1/2 restriction. 
  #' 
  #' @param in_df A dataframe-output of df_obs_logbase_quartiles_info() or df_obs_quartiles_info().
  #' @param N_obs_bin Threshold on number of observations per bin.
  #' @param morethan Theshold on gene coverage for lm (default = 10).
  #' @return log2(lm intercept).
  #' @examples
  #' 
  df = in_df[in_df$MeanGeneCoverage > morethan & in_df$N_observations > N_obs_bin, 
             c('MeanGeneCoverage','deltaAI','N_observations')]
  df[, c('MeanGeneCoverage','deltaAI')] = log2(df[, c('MeanGeneCoverage','deltaAI')])
  lm(data = df, deltaAI ~ offset(-0.5*MeanGeneCoverage))$coefficients[1]
}

# .......................................................................................
# Aftercomments:
# _______________________________________________________________________________________


# THE END