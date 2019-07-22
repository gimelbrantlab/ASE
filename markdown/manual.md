Allele-specific analysis of RNA-seq data using technical replicates
================

We start by loading the libraries and functions we will need:

```r
library("tidyverse")
library("knitr")
library("biomaRt")
source("R/ASE_functions.R")
source("R/PerformAIAnalysis_CC.R")
```

### Data

To demonstrate how the pipeline works, we are going to use data from RNA-seq experiment aimed to study effect of 5aza treatment on expression and allele-specific expression. Let's look at the data first:

### Experiment design

Abelson clones:

| #       | clone                         | number of repicates |
|---------|-------------------------------|---------------------|
| 1       | 1.1                           | 2                   |
| 2       | 4.11                          | 2                   |
| 3       | H8                            | 2                   |

The first step is to create design matrix describing the experiment:

``` r
experimentNames <- c("clone_1.1","clone_4.11","clone_H8")
techReps <- c(2,2,2)
designMatrix <- BuildDesign(experimentNames, techReps)
```

### Load the data

Now we can look at the data. Here is the *geneCountTab* dataframe with counts for all genes, all replicates and conditions:

``` r
inTabs <- "~/Dropbox (Partners HealthCare)/replicates_ASE/data/5aza/5aza_no_treatment_3clones_processed_gene_extended2.txt"
geneCountTab <- GetGatkPipelineTabs(inTabs, c(2,2,2), multiple = F, chrom = F)
colnames(geneCountTab)[2:13] <- c(nameColumns(1,2), nameColumns(2,2), nameColumns(3,2))
head(geneCountTab)
```


### Look at AI correlations

We can calculate allelic imbalances for all genes for all experimental conditions pooling technical replicates:
$$AI=\\frac{maternal\\ counts}{maternal\\ counts + paternal\\ counts}$$

``` r
aiTable <- do.call(cbind, lapply(1:length(designMatrix$techReps), function(x){
  round(CountsToAI(geneCountTab, reps = unlist(designMatrix$replicateNums[x]),thr=10)$AI,3)
}))
aiTable <- data.frame(geneCountTab$ensembl_gene_id, aiTable)
colnames(aiTable) <- c("ensembl_gene_id", designMatrix$experimentNames)
head(na.omit(aiTable))
```

    ## ensembl_gene_id clone_1.1 clone_4.11 clone_H8
    ## 8  ENSMUSG00000033845     0.575      0.530    0.604
    ## 9  ENSMUSG00000025903     0.531      0.535    0.524
    ## 10 ENSMUSG00000033813     0.555      0.555    0.565
    ## 11 ENSMUSG00000062588     0.010      0.026    0.007
    ## 15 ENSMUSG00000033793     0.458      0.518    0.467
    ## 18 ENSMUSG00000090031     0.442      0.387    0.163

Let's visualize AI correlations:

``` r
ggplot(aiTable, aes_string(x=as.name(designMatrix$experimentNames[1]), y=as.name(designMatrix$experimentNames[2]))) +
  geom_point(size=0.5) +
  theme_bw() +
  coord_fixed() 
```


### Differential AI analysis

Let's compare AIs in two Abelson clones . We will construct 95%-CIs around AIs and get the resulting classification finding the genes demonstrating difference in AI between clones (based on non-overlapping CIs).

First, we need to get QCC for the experiments:



``` r
CIs_all <- lapply(1:length(designMatrix$techReps), function(x){
   resCC <- PerformBinTestAIAnalysisForConditionNPoint(geneCountTab,
                                                       unlist(designMatrix$replicateNums[x]),
                                                       Q=0.95,
                                                       EPS=1.3,
                                                       thr=NA)
   resCC$CC
 })
sink("~/Dropbox (Partners HealthCare)/replicates_ASE/data/5aza/5aza_no_treatment_3clones_CC_backup_upd2.txt")
writeLines(unlist(lapply(CIs_all, paste, collapse=" ")))
sink()
```
Alternatively, if we've ran this analysis before, we can load pre-computer QCC estimates:

``` r
CC_backup <- read.table("~/Dropbox (Partners HealthCare)/replicates_ASE/data/5aza/5aza_no_treatment_3clones_CC_backup_upd2.txt", header=FALSE, fill = TRUE)
CC_backup
```

Now let's run the differential analysis. We define a gene as demonstrating differential allelic imbalance between two clones or conditions if this gene has non-intersecting CIs of allelic imbalance estimates in these clones/conditions and if the difference between these allelic imbalance estimates is bigger than some cutoff (set to 0.1 here). The cutoff value is set arbitary depending on the biological question we are asking.


``` r
thr_coverage <- 40
minDifference <- 0.1
experimentA <- 1
experimentB <- 2

df_compare <- PerformBinTestAIAnalysisForTwoConditions_knownCC(geneCountTab, 
                                                               vect1CondReps = unlist(designMatrix$replicateNums[experimentA]), 
                                                               vect2CondReps = unlist(designMatrix$replicateNums[experimentB]), 
                                                               vect1CondRepsCombsCC = as.numeric(CC_backup[experimentA,!is.na(CC_backup[experimentA,])]),
                                                               vect2CondRepsCombsCC = as.numeric(CC_backup[experimentB,!is.na(CC_backup[experimentB,])]),
                                                               Q = 0.95,
                                                               thr = thr_coverage,
                                                               minDifference = minDifference
)
AI_table_with_CIs <- df_compare[,c(1,4,21,22,11,23,24,25,26)]
colnames(AI_table_with_CIs) <- c("ensembl_gene_id", 
                                 paste0("AI_",designMatrix$experimentNames[experimentA]),
                                 paste0("AI_CI_left_",designMatrix$experimentNames[experimentA]), 
                                 paste0("AI_CI_right_",designMatrix$experimentNames[experimentA]), 
                                 paste0("AI_",designMatrix$experimentNames[experimentB]),
                                 paste0("AI_CI_left_",designMatrix$experimentNames[experimentB]), 
                                 paste0("AI_CI_right_",designMatrix$experimentNames[experimentB]), 
                                 "CI_diff", 
                                 "CI_plus_minDiff")
```

We can plot resulting analysis and color genes based on differential allelic expression:

``` r
fig_compare <- ggplot(df_compare, aes(x = AI_1, y = AI_2, col = BT_CC_thrDiff)) + 
  geom_point(size=0.5) + 
  theme_bw() +
  xlab(paste0("AI, ",designMatrix$experimentNames[experimentA])) +
  ylab(paste0("AI, ",designMatrix$experimentNames[experimentB])) +
  scale_color_manual(name="Differential AI", labels=c("FALSE", "TRUE"), values=c("gray", "red")) +
  coord_fixed()
fig_compare
```

Next, if we want to look at numbers, it is useful to add chromosome information:

``` r
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
genes <- getBM(c("ensembl_gene_id", "chromosome_name"), mart=ensembl)
genes$chr <- paste0("chr", genes$chromosome_name)
AI_table_with_CIs <- merge(AI_table_with_CIs, genes[,c("ensembl_gene_id", "chr")], all.x = "T", by = "ensembl_gene_id")
``` 

``` r
print(paste0("Number of genes demonstrating differential allelic imbalance: ", table(AI_table_with_CIs$CI_plus_minDiff)["TRUE"], " out of ", sum(!is.na(AI_table_with_CIs$CI_plus_minDiff))))
```

    ## [1] "Number of genes demonstrating differential allelic imbalance: 815 out of 9121"

Now if we exclude X chromosome:

``` r
print(paste0("Number of genes demonstrating differential allelic imbalance: ", table(AI_table_with_CIs$CI_plus_minDiff[AI_table_with_CIs$chr!="chrX" & AI_table_with_CIs$chr!="chrY"])["TRUE"], " out of ", sum(!is.na(AI_table_with_CIs$CI_plus_minDiff[AI_table_with_CIs$chr!="chrX" & AI_table_with_CIs$chr!="chrY"]))))
```

    ## [1] "Number of autosomal genes demonstrating differential allelic imbalance: 606 out of 8815"

