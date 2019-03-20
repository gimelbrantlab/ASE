# load data (sampled 50 mln reads, processed)
tabfile <- "~/Dropbox (Partners HealthCare)/replicates_ASE/data/kidney/submln/SMART10ng_mln50_processed_gene_extended2.txt"
df <- GetGatkPipelineTabs(tabfile, nReps = c(30), multiple=F, chrom = T)
nameColumns <- function(exp_n, rep_n)  {
  paste0("exp", rep(exp_n, 2*rep_n), "_rep", rep(1:rep_n, each = 2), "_", rep(c("ref", "alt"), rep_n))
}
colnames(df)[2:(length(colnames(df))-1)] <- c(nameColumns(1,5), nameColumns(2,5), nameColumns(3,5), nameColumns(4,5), nameColumns(5,5), nameColumns(6,5))
df <- df %>% filter(chr != "chrX") %>% dplyr::select(-chr)
# sample experiments
samplereps <- (0:5)*5 + sample(1:5, 6, rep=T)
df_sample <- df[, sort(c(1, samplereps*2, samplereps*2+1))]
colnames(df_sample)[(1:6)*2] <- paste0(rep("rep_",6), 1:6, "_ref")
colnames(df_sample)[1+(1:6)*2] <- paste0(rep("rep_",6), 1:6, "_alt")
df_sample <- df_sample %>% 
  mutate(rep_all_ref = rep_1_ref + rep_2_ref + rep_3_ref + rep_4_ref + rep_5_ref + rep_6_ref) %>% 
  mutate(rep_all_alt = rep_1_alt + rep_2_alt + rep_3_alt + rep_4_alt + rep_5_alt + rep_6_alt) 
# get AIs
aiTable <- do.call(cbind, lapply(1:6, function(x){
  round(CountsToAI(df_sample, reps = x, thr=10)$AI,3)
}))
aiTable6 <- round(CountsToAI(df_sample, reps = 1:6, thr=10)$AI,3)
aiTable12 <- round(CountsToAI(df_sample, reps = 1:2, thr=10)$AI,3)
aiTable34 <- round(CountsToAI(df_sample, reps = 3:4, thr=10)$AI,3)
aiTable <- data.frame(df_sample$ensembl_gene_id, aiTable, aiTable6, aiTable12, aiTable34)
colnames(aiTable) <- c("ensembl_gene_id", paste0(rep("rep_",6), 1:6, "AI"), "all6_AI", "rep_12AI", "rep_34AI")
head(aiTable)
#
# get coverages
covTable <- do.call(cbind, lapply(1:6, function(x){
  MeanCoverage(df_sample, reps = x)$meanCOV
}))
covTable6 <- round(MeanCoverage(df_sample, reps = 1:6)$meanCOV)
covTable <- data.frame(df_sample$ensembl_gene_id, covTable, covTable6)
colnames(covTable) <- c("ensembl_gene_id", paste0(rep("rep_",6), 1:6, "cov"), "all6_cov")
head(covTable)
# merge tables
finalTable <- merge(aiTable, covTable, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
# get binomial pvals
binom.test.p <- function(x) {
  if (is.na(x[1])|is.na(x[2])|(x[1]+x[2])<10) {
    return(NA)
  } 
  else {
    return(binom.test(x, alternative="two.sided",conf.level = 0.95, p=0.5)$p.value)
  }
} 
df_sample_binomial_pvals <- df_sample %>% mutate(rep1_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_1_ref, rep_1_alt)) %>% 
  mutate(rep2_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_2_ref, rep_2_alt)) %>% 
  mutate(rep3_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_3_ref, rep_3_alt)) %>% 
  mutate(rep4_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_4_ref, rep_4_alt)) %>% 
  mutate(rep5_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_5_ref, rep_5_alt)) %>% 
  mutate(rep6_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_6_ref, rep_6_alt)) %>% 
  mutate(rep_all_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_all_ref, rep_all_alt)) %>% 
  mutate(rep1_BT0.95 = rep1_pval < 0.05/nrow(finalTable)) %>% 
  mutate(rep2_BT0.95 = rep2_pval < 0.05/nrow(finalTable)) %>% 
  mutate(rep3_BT0.95 = rep3_pval < 0.05/nrow(finalTable)) %>% 
  mutate(rep4_BT0.95 = rep4_pval < 0.05/nrow(finalTable)) %>% 
  mutate(rep5_BT0.95 = rep5_pval < 0.05/nrow(finalTable)) %>% 
  mutate(rep6_BT0.95 = rep6_pval < 0.05/nrow(finalTable)) %>% 
  mutate(all6_BT0.95 = rep_all_pval < 0.05/nrow(finalTable))
df_sample_binomial_pvals <- df_sample_binomial_pvals %>% mutate(rep12_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_1_ref+rep_2_ref, rep_1_alt+rep_2_alt)) %>% 
  mutate(rep34_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_3_ref+rep_4_ref, rep_3_alt+rep_4_alt)) %>% 
  mutate(rep12_BT0.95 = rep12_pval < 0.05/nrow(finalTable)) %>% 
  mutate(rep34_BT0.95 = rep34_pval < 0.05/nrow(finalTable))
# can out method results
rep12_CIs <- PerformDiffAIAnalysisForConditionNPoint(df_sample, 
                                                          vectReps = 1:2,
                                                          condName="Condition", 
                                                          pt = 0.5,
                                                          thr = 40,
                                                     minDifference = 0.1)
rep34_CIs <- PerformDiffAIAnalysisForConditionNPoint(df_sample, 
                                                     vectReps = 3:4,
                                                     condName="Condition", 
                                                     pt = 0.5,
                                                     thr = 40,
                                                     minDifference = 0.1)
CI_method_res <- merge(rep12_CIs[,c(1,3,4,9)], rep34_CIs[,c(1,3,4,9)], by.x="ID", by.y="ID")
colnames(CI_method_res) <- c("ensembl_gene_id", "meanCov_12", "meanAI_12", "rep12_diffAI", "meanCov_34", "meanAI_34","rep34_diffAI")
# merge tables
finalTable_plus <- merge(finalTable, df_sample_binomial_pvals[,c(1,16:33)], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
finalTable_plus <- merge(finalTable_plus, CI_method_res, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
finalTable_plus <- na.omit(finalTable_plus)
# plots
grid.arrange(
  ggplot(finalTable_plus, aes(rep_1cov, rep_1AI, color = rep1_BT0.95)) +
    geom_point(size=0.2) +
    theme_linedraw() +
    facet_grid(~rep1_BT0.95) +
    scale_x_continuous(trans='log2'),
  ggplot(finalTable_plus, aes(rep_2cov, rep_2AI, color = rep2_BT0.95)) +
    geom_point(size=0.2) +
    theme_linedraw() +
    facet_grid(~rep1_BT0.95) +
    scale_x_continuous(trans='log2'),
  ggplot(finalTable_plus, aes(all6_cov, all6_AI, color = all6_BT0.95)) +
    geom_point(size=0.2) +
    theme_linedraw() +
    facet_grid(~rep1_BT0.95) +
    scale_x_continuous(trans='log2')
)
#
grid.arrange(
  ggplot(finalTable_plus, aes(rep_1cov+rep_2cov, rep_12AI, color = rep12_BT0.95)) +
    geom_point(size=0.2) +
    theme_linedraw() +
    facet_grid(~rep12_BT0.95) +
    scale_x_continuous(trans='log2'),
  ggplot(finalTable_plus, aes(rep_3cov+rep_4cov, rep_34AI, color = rep34_BT0.95)) +
    geom_point(size=0.2) +
    theme_linedraw() +
    facet_grid(~rep12_BT0.95) +
    scale_x_continuous(trans='log2')
)
#
grid.arrange(
  ggplot(CI_method_res, aes(meanCov_12, meanAI_12, color = rep12_diffAI)) +
    geom_point(size=0.2) +
    theme_linedraw() +
    facet_grid(~rep12_diffAI) +
    scale_x_continuous(trans='log2'),
  ggplot(CI_method_res, aes(meanCov_34, meanAI_34, color = rep34_diffAI)) +
    geom_point(size=0.2) +
    theme_linedraw() +
    facet_grid(~rep12_diffAI) +
    scale_x_continuous(trans='log2')
)


