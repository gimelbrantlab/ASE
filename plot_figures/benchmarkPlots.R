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
aiTable <- data.frame(df_sample$ensembl_gene_id, aiTable, aiTable6)
colnames(aiTable) <- c("ensembl_gene_id", paste0(rep("rep_",6), 1:6, "AI"), "all6_AI")
head(aiTable)
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
df_sample_pvals <- df_sample %>% mutate(rep1_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_1_ref, rep_1_alt)) %>% 
  mutate(rep2_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_2_ref, rep_2_alt)) %>% 
  mutate(rep3_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_3_ref, rep_3_alt)) %>% 
  mutate(rep4_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_4_ref, rep_4_alt)) %>% 
  mutate(rep5_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_5_ref, rep_5_alt)) %>% 
  mutate(rep6_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_6_ref, rep_6_alt)) %>% 
  mutate(rep_all_pval = mapply(function(x, y) binom.test.p(c(x,y)), rep_all_ref, rep_all_alt))
# merge tables
finalTable <- merge(finalTable, df_sample_pvals[,c(1,16:22)], by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
# add binomial test outputs
finalTable <- finalTable %>% mutate(rep1_BT0.95 = rep1_pval < 0.05/nrow(finalTable)) %>% 
  mutate(rep2_BT0.95 = rep2_pval < 0.05/nrow(finalTable)) %>% 
  mutate(rep3_BT0.95 = rep3_pval < 0.05/nrow(finalTable)) %>% 
  mutate(rep4_BT0.95 = rep4_pval < 0.05/nrow(finalTable)) %>% 
  mutate(rep5_BT0.95 = rep5_pval < 0.05/nrow(finalTable)) %>% 
  mutate(rep6_BT0.95 = rep6_pval < 0.05/nrow(finalTable)) %>% 
  mutate(all6_BT0.95 = rep_all_pval < 0.05/nrow(finalTable))
finalTable <- na.omit(finalTable)
# plots
grid.arrange(
  ggplot(finalTable, aes(rep_1cov, rep_1AI, color = rep1_BT0.95)) +
    geom_point(size=0.2) +
    theme_linedraw() +
    facet_grid(~rep1_BT0.95) +
    scale_x_continuous(trans='log2'),
  ggplot(finalTable, aes(rep_2cov, rep_2AI, color = rep2_BT0.95)) +
    geom_point(size=0.2) +
    theme_linedraw() +
    facet_grid(~rep1_BT0.95) +
    scale_x_continuous(trans='log2'),
  ggplot(finalTable, aes(all6_cov, all6_AI, color = all6_BT0.95)) +
    geom_point(size=0.2) +
    theme_linedraw() +
    facet_grid(~rep1_BT0.95) +
    scale_x_continuous(trans='log2')
)
