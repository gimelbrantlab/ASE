setwd("/home/asya/Dropbox (Partners HealthCare)/code/ASE/plot_figures")

library(tidyverse)
library(fitdistrplus)
library(cowplot)
library(ggrepel)
library(VennDiagram)

source("../R/ASE_functions.R")
source("../R/PerformAIAnalysis_CC.R")
source("../R/DownstreamASE.R")
source("../plot_figures/utilPlots.R")

MLN = 30

# Load X-mln data:
NEB_data = GetGatkPipelineTabs(paste0("../../../data/kidney/submln/", "MLN", MLN,"_SG", 1:6, "_N955_", "NEB", "_R1_merged_v2.mln", MLN,"_trial5_processed_gene_extended2.txt"), c(5,5,5,5,5,5), multiple = T)
SMART10ng_data = GetGatkPipelineTabs(paste0("../../../data/kidney/submln/", "MLN", MLN,"_SG", 7:12, "_N955_", "SMARTseq_10_ng", "_R1_merged_v2.mln", MLN,"_trial5_processed_gene_extended2.txt"), c(5,5,5,5,5,5), multiple = T)
SMART100pg_data = GetGatkPipelineTabs(paste0("../../../data/kidney/submln/", "MLN", MLN,"_SG", 1:6, "_N955_", "SMARTseq_100_pg", "_R1_merged_v2.mln", MLN,"_trial5_processed_gene_extended2.txt"), c(5,5,5,5,5,5), multiple = T)

data = list(NEB_data, SMART10ng_data, SMART100pg_data)
rm(NEB_data)
rm(SMART10ng_data)
rm(SMART100pg_data)

removeX <- function(DF, legitim_chrgenes){
  return(DF[DF$ensembl_gene_id %in% legitim_chrgenes$gene, ])
}
chrgenes <- read.delim("../../../data/kidney/submln/MLN30_SG3_N955_NEB_R1_merged_v2.mln30_trial5_processed_gene_extended2.txt")[, c("ensembl_gene_id", "chr")]

data_noX = lapply(data, function(x){
  x[x$ensembl_gene_id %in% chrgenes[chrgenes$chr!="chrX" & chrgenes$chr!="chrY", "ensembl_gene_id"], ]
})

CC = lapply(c(paste0("../../../data/kidney/submln/NEB_CorrConsts_", MLN,"mln_1.05.RData"),
              paste0("../../../data/kidney/submln/SMARTseq_10_ng_CorrConsts_", MLN,"mln_1.05.RData"),
              paste0("../../../data/kidney/submln/SMARTseq_100_pg_CorrConsts_", MLN,"mln_1.05.RData")),
            function(file){
              load(file)
              sapply(out_XXmln_SMART10ng, function(x){x$fittedCC})
            })

# DFs for Euler diagrams:

CalculateTestWithHalf <- function(df, cc, exp_name, thr=10){
  res = Reduce(function(x, y) merge(x, y, by="ID"), 
               lapply(1:(ncol(df)%/%(2*5)), function(j){
                 sample1rep = (j-1)*5 + sample(1:5, 1)
                 PerformBinTestAIAnalysisForConditionNPoint_knownCC(df, sample1rep, cc, thr=thr)[, c("ID", "BT", "BT_CC")]
               })
  )
  names(res) = c("ID", 
                 paste0(c("BT.", "BT_CC."), 
                        rep(1:(ncol(df)%/%(2*5)), each=2)
                 )
  )
  return(res)
}


# PLOTTING:
# (TRUE -- biased, FALSE -- not)
# BT.x is res of Binomial test (Bonf corr) for x-th replicate
# BT_CC.x is res of Corrected on CC Binomial test for x-th replicate

# example SMART 10 ng:
df = CalculateTestWithHalf(data_noX[[2]], list_of_consts[[2]], list_of_libprepnames[[2]], thr=10)

# Pair:
# example for 1,2,3: Binomial and CC corrected Binomial
grid.newpage()
draw.pairwise.venn(nrow(subset(df, BT_CC.1 == "TRUE")), 
                   nrow(subset(df, BT_CC.2 == "TRUE")), 
                   nrow(subset(df, BT_CC.1 == "TRUE" & BT_CC.2 == "TRUE")), 
                   category = c("Tech Replicate #1", "Tech Replicate #2"), 
                   lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
                   scaled = TRUE, print.mode=c("raw", "percent"))
grid.newpage()
draw.pairwise.venn(nrow(subset(df, BT.1 == "TRUE")), 
                   nrow(subset(df, BT.2 == "TRUE")), 
                   nrow(subset(df, BT.1 == "TRUE" & BT_CC.2 == "TRUE")), 
                   category = c("Tech Replicate #1", "Tech Replicate #2"), 
                   lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
                   scaled = TRUE, print.mode=c("raw", "percent"))
# Triple:
# example for 1,2,3
grid.newpage()
draw.triple.venn(nrow(subset(df, BT.1 == "TRUE")), 
                 nrow(subset(df, BT.2 == "TRUE")), 
                 nrow(subset(df, BT.3 == "TRUE")),
                 nrow(subset(df, BT.1 == "TRUE" & BT.2 == "TRUE")), 
                 nrow(subset(df, BT.1 == "TRUE" & BT.3 == "TRUE")),
                 nrow(subset(df, BT.2 == "TRUE" & BT.3 == "TRUE")),
                 nrow(subset(df, BT.1 == "TRUE" & BT.2 == "TRUE" & BT.3 == "TRUE")),
                 category = c("Tech Replicate #1", "Tech Replicate #2", "Tech Replicate #3"), 
                 lty = rep("blank", 3), fill = c("light blue", "pink", "light yellow"), alpha = rep(0.5, 3), 
                 scaled = TRUE, print.mode=c("raw", "percent"))
grid.newpage()
draw.triple.venn(nrow(subset(df, BT_CC.1 == "TRUE")), 
                 nrow(subset(df, BT_CC.2 == "TRUE")), 
                 nrow(subset(df, BT_CC.3 == "TRUE")),
                 nrow(subset(df, BT_CC.1 == "TRUE" & BT_CC.2 == "TRUE")), 
                 nrow(subset(df, BT_CC.1 == "TRUE" & BT_CC.3 == "TRUE")),
                 nrow(subset(df, BT_CC.2 == "TRUE" & BT_CC.3 == "TRUE")),
                 nrow(subset(df, BT_CC.1 == "TRUE" & BT_CC.2 == "TRUE" & BT_CC.3 == "TRUE")),
                 category = c("Tech Replicate #1", "Tech Replicate #2", "Tech Replicate #3"), 
                 lty = rep("blank", 3), fill = c("light blue", "pink", "light yellow"), alpha = rep(0.5, 3), 
                 scaled = TRUE, print.mode=c("raw", "percent"))
# Quadro:
# example for 1,2,3,4
grid.newpage()
draw.quad.venn(nrow(subset(df, BT_CC.1 == "TRUE")), 
               nrow(subset(df, BT_CC.2 == "TRUE")), 
               nrow(subset(df, BT_CC.3 == "TRUE")),
               nrow(subset(df, BT_CC.4 == "TRUE")),
               nrow(subset(df, BT_CC.1 == "TRUE" & BT_CC.2 == "TRUE")), 
               nrow(subset(df, BT_CC.1 == "TRUE" & BT_CC.3 == "TRUE")),
               nrow(subset(df, BT_CC.1 == "TRUE" & BT_CC.4 == "TRUE")),
               nrow(subset(df, BT_CC.2 == "TRUE" & BT_CC.3 == "TRUE")),
               nrow(subset(df, BT_CC.2 == "TRUE" & BT_CC.4 == "TRUE")),
               nrow(subset(df, BT_CC.3 == "TRUE" & BT_CC.4 == "TRUE")),
               nrow(subset(df, BT_CC.1 == "TRUE" & BT_CC.2 == "TRUE" & BT_CC.3 == "TRUE")),
               nrow(subset(df, BT_CC.1 == "TRUE" & BT_CC.2 == "TRUE" & BT_CC.4 == "TRUE")),
               nrow(subset(df, BT_CC.1 == "TRUE" & BT_CC.3 == "TRUE" & BT_CC.4 == "TRUE")),
               nrow(subset(df, BT_CC.2 == "TRUE" & BT_CC.3 == "TRUE" & BT_CC.4 == "TRUE")),
               nrow(subset(df, BT_CC.1 == "TRUE" & BT_CC.2 == "TRUE" & BT_CC.3 == "TRUE" & BT_CC.4 == "TRUE")),
               category = c("Tech Replicate #1", "Tech Replicate #2", "Tech Replicate #3", "Tech Replicate #4"), 
               lty = rep("blank", 4), fill = c("light blue", "pink", "light yellow", "light green"), alpha = rep(0.5, 4), 
               scaled = TRUE, print.mode=c("raw", "percent"))
# Quadro:
# example for 1,2 and comparing Binomial with CC corrected Binomial
grid.newpage()
draw.quad.venn(nrow(subset(df, BT.1 == "TRUE")), 
               nrow(subset(df, BT.2 == "TRUE")), 
               nrow(subset(df, BT_CC.1 == "TRUE")),
               nrow(subset(df, BT_CC.2 == "TRUE")),
               nrow(subset(df, BT.1 == "TRUE" & BT.2 == "TRUE")), 
               nrow(subset(df, BT.1 == "TRUE" & BT_CC.1 == "TRUE")),
               nrow(subset(df, BT.1 == "TRUE" & BT_CC.2 == "TRUE")),
               nrow(subset(df, BT.2 == "TRUE" & BT_CC.1 == "TRUE")),
               nrow(subset(df, BT.2 == "TRUE" & BT_CC.2 == "TRUE")),
               nrow(subset(df, BT_CC.1 == "TRUE" & BT_CC.2 == "TRUE")),
               nrow(subset(df, BT.1 == "TRUE" & BT.2 == "TRUE" & BT_CC.1 == "TRUE")),
               nrow(subset(df, BT.1 == "TRUE" & BT.2 == "TRUE" & BT_CC.2 == "TRUE")),
               nrow(subset(df, BT.1 == "TRUE" & BT_CC.1 == "TRUE" & BT_CC.2 == "TRUE")),
               nrow(subset(df, BT.2 == "TRUE" & BT_CC.1 == "TRUE" & BT_CC.2 == "TRUE")),
               nrow(subset(df, BT.1 == "TRUE" & BT.2 == "TRUE" & BT_CC.1 == "TRUE" & BT_CC.2 == "TRUE")),
               category = c("Tech Replicate #1\n(BT)", "Tech Replicate #2\n(BT)", "Tech Replicate #1\n(BT_CC)", "Tech Replicate #2\n(BT_CC)"), 
               lty = rep("blank", 4), fill = c("light blue", "pink", "light yellow", "light green"), alpha = rep(0.5, 4), 
               scaled = TRUE, print.mode=c("raw", "percent"))


# Less handly:

df_NEB = CalculateTestWithHalf(data_noX[[1]], list_of_consts[[1]], list_of_libprepnames[[1]], thr=10)
df_SM10 = CalculateTestWithHalf(data_noX[[2]], list_of_consts[[2]], list_of_libprepnames[[2]], thr=10)
df_SM100 = CalculateTestWithHalf(data_noX[[3]], list_of_consts[[3]], list_of_libprepnames[[3]], thr=10)

df_NEB_SM10_SM100 = Reduce(function(x, y) merge(x, y, by="ID"), list(df_NEB, df_SM10, df_SM100))


grid.newpage(); draw.pairwise.venn(nrow(subset(df_NEB_SM10_SM100, BT.2.x == "TRUE")), 
                   nrow(subset(df_NEB_SM10_SM100, BT.2.y == "TRUE")), 
                   nrow(subset(df_NEB_SM10_SM100, BT.2.x == "TRUE" & BT.2.y == "TRUE")), 
                   category = c("Tech Replicate \nNEB", "Tech Replicate \nSMART10ng"), 
                   lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
                   scaled = TRUE, print.mode=c("raw", "percent"))
grid.newpage(); draw.pairwise.venn(nrow(subset(df_NEB_SM10_SM100, BT.2.x == "TRUE")), 
                   nrow(subset(df_NEB_SM10_SM100, BT.2 == "TRUE")), 
                   nrow(subset(df_NEB_SM10_SM100, BT.2.x == "TRUE" & BT.2 == "TRUE")), 
                   category = c("Tech Replicate \nNEB", "Tech Replicate \nSMART100pg"), 
                   lty = rep("blank", 2), fill = c("light blue", "light green"), alpha = rep(0.5, 2), 
                   scaled = TRUE, print.mode=c("raw", "percent"))
grid.newpage(); draw.pairwise.venn(nrow(subset(df_NEB_SM10_SM100, BT.2.y == "TRUE")), 
                   nrow(subset(df_NEB_SM10_SM100, BT.2 == "TRUE")), 
                   nrow(subset(df_NEB_SM10_SM100, BT.2.y == "TRUE" & BT.2 == "TRUE")), 
                   category = c("Tech Replicate \nSMART10ng", "Tech Replicate \nSMART100pg"), 
                   lty = rep("blank", 2), fill = c("pink", "light green"), alpha = rep(0.5, 2), 
                   scaled = TRUE, print.mode=c("raw", "percent"))


grid.newpage(); draw.pairwise.venn(nrow(subset(df_NEB_SM10_SM100, BT.2.x == "TRUE")), 
                   nrow(subset(df_NEB_SM10_SM100, BT.3.x == "TRUE")), 
                   nrow(subset(df_NEB_SM10_SM100, BT.2.x == "TRUE" & BT.3.x == "TRUE")), 
                   category = c("Tech Replicate \nNEB", "Tech Replicate \nNEB"), 
                   fill = c("light blue", "light blue"), alpha = rep(0.5, 2), 
                   scaled = TRUE, print.mode=c("raw", "percent"))
grid.newpage(); draw.pairwise.venn(nrow(subset(df_NEB_SM10_SM100, BT.2.y == "TRUE")), 
                   nrow(subset(df_NEB_SM10_SM100, BT.3.y == "TRUE")), 
                   nrow(subset(df_NEB_SM10_SM100, BT.2.y == "TRUE" & BT.3.y == "TRUE")), 
                   category = c("Tech Replicate \nSMART10ng", "Tech Replicate \nSMART10ng"), 
                   fill = c("pink", "pink"), alpha = rep(0.5, 2), 
                   scaled = TRUE, print.mode=c("raw", "percent"))
grid.newpage(); draw.pairwise.venn(nrow(subset(df_NEB_SM10_SM100, BT.2 == "TRUE")), 
                   nrow(subset(df_NEB_SM10_SM100, BT.3 == "TRUE")), 
                   nrow(subset(df_NEB_SM10_SM100, BT.2 == "TRUE" & BT.3 == "TRUE")), 
                   category = c("Tech Replicate \nSMART100pg", "Tech Replicate \nSMART100pg"), 
                   fill = c("light green", "light green"), alpha = rep(0.5, 2), 
                   scaled = TRUE, print.mode=c("raw", "percent"))


grid.newpage(); draw.pairwise.venn(nrow(subset(df_NEB_SM10_SM100, BT_CC.2.x == "TRUE")), 
                                   nrow(subset(df_NEB_SM10_SM100, BT_CC.2.y == "TRUE")), 
                                   nrow(subset(df_NEB_SM10_SM100, BT_CC.2.x == "TRUE" & BT_CC.2.y == "TRUE")), 
                                   category = c("Tech Replicate \nNEB", "Tech Replicate \nSMART10ng"), 
                                   lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
                                   scaled = TRUE, print.mode=c("raw", "percent"))
grid.newpage(); draw.pairwise.venn(nrow(subset(df_NEB_SM10_SM100, BT_CC.2.x == "TRUE")), 
                                   nrow(subset(df_NEB_SM10_SM100, BT_CC.2 == "TRUE")), 
                                   nrow(subset(df_NEB_SM10_SM100, BT_CC.2.x == "TRUE" & BT_CC.2 == "TRUE")), 
                                   category = c("Tech Replicate \nNEB", "Tech Replicate \nSMART100pg"), 
                                   lty = rep("blank", 2), fill = c("light blue", "light green"), alpha = rep(0.5, 2), 
                                   scaled = TRUE, print.mode=c("raw", "percent"))
grid.newpage(); draw.pairwise.venn(nrow(subset(df_NEB_SM10_SM100, BT_CC.2.y == "TRUE")), 
                                   nrow(subset(df_NEB_SM10_SM100, BT_CC.2 == "TRUE")), 
                                   nrow(subset(df_NEB_SM10_SM100, BT_CC.2.y == "TRUE" & BT_CC.2 == "TRUE")), 
                                   category = c("Tech Replicate \nSMART10ng", "Tech Replicate \nSMART100pg"), 
                                   lty = rep("blank", 2), fill = c("pink", "light green"), alpha = rep(0.5, 2), 
                                   scaled = TRUE, print.mode=c("raw", "percent"))


grid.newpage(); draw.pairwise.venn(nrow(subset(df_NEB_SM10_SM100, BT_CC.2.x == "TRUE")), 
                                   nrow(subset(df_NEB_SM10_SM100, BT_CC.3.x == "TRUE")), 
                                   nrow(subset(df_NEB_SM10_SM100, BT_CC.2.x == "TRUE" & BT_CC.3.x == "TRUE")), 
                                   category = c("Tech Replicate \nNEB", "Tech Replicate \nNEB"), 
                                   fill = c("light blue", "light blue"), alpha = rep(0.5, 2), 
                                   scaled = TRUE, print.mode=c("raw", "percent"))
grid.newpage(); draw.pairwise.venn(nrow(subset(df_NEB_SM10_SM100, BT_CC.2.y == "TRUE")), 
                                   nrow(subset(df_NEB_SM10_SM100, BT_CC.3.y == "TRUE")), 
                                   nrow(subset(df_NEB_SM10_SM100, BT_CC.2.y == "TRUE" & BT_CC.3.y == "TRUE")), 
                                   category = c("Tech Replicate \nSMART10ng", "Tech Replicate \nSMART10ng"), 
                                   fill = c("pink", "pink"), alpha = rep(0.5, 2), 
                                   scaled = TRUE, print.mode=c("raw", "percent"))
grid.newpage(); draw.pairwise.venn(nrow(subset(df_NEB_SM10_SM100, BT_CC.2 == "TRUE")), 
                                   nrow(subset(df_NEB_SM10_SM100, BT_CC.3 == "TRUE")), 
                                   nrow(subset(df_NEB_SM10_SM100, BT_CC.2 == "TRUE" & BT_CC.3 == "TRUE")), 
                                   category = c("Tech Replicate \nSMART100pg", "Tech Replicate \nSMART100pg"), 
                                   fill = c("light green", "light green"), alpha = rep(0.5, 2), 
                                   scaled = TRUE, print.mode=c("raw", "percent"))






# # Subsample 2 reps from data:
# 
# Ntrials = 1
# set.seed(1)
# sampleMreps = sapply(1:(max(Ntrials, 10)), function(i){(sample(1:6, 2)-1)*5 + sample(5, 2, replace=T)})
# list_of_datas = data_noX
# list_of_consts = list(mean(CC[[1]][6:15]), mean(CC[[2]]), mean(CC[[3]]))
# list_of_libprepnames = list("1: NEBNext (100ng)", "2: SMARTseq (10ng)", "3: SMARTseq (0.1ng)")
# reppair = sampleMreps[, 1]
# 
# DF_forplot = CreateForplotDF_btNbtcc(lapply(1:3, function(i){
#   CreateForplotDF(list_of_datas[[i]], reppair, list_of_consts[[i]], list_of_libprepnames[[i]])
# }))
# 
# DF_forplot$BTBFCC$BT.x = as.factor(DF_forplot$BTBFCC$BT.x)
# levels(DF_forplot$BTBFCC$BT.x) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")
# DF_forplot$BTBFCC$BT.y = as.factor(DF_forplot$BTBFCC$BT.y)
# levels(DF_forplot$BTBFCC$BT.y) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")
# 
# DF_forplot$BTBF$BT.x = as.factor(DF_forplot$BTBF$BT.x)
# levels(DF_forplot$BTBF$BT.x) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")
# DF_forplot$BTBF$BT.y = as.factor(DF_forplot$BTBF$BT.y)
# levels(DF_forplot$BTBF$BT.y) = c("Balanced genes \n (according to replicate 1)", "Imbalanced Genes \n (according to replicate 1)")
# 
# # Venn-Euler for given pairs:
# # No CC
# for(lib in unique(DF_forplot$BTBF$libprep)){
#   df = DF_forplot$BTBF[DF_forplot$BTBF$libprep == lib, ]
#   grid.newpage()
#   draw.pairwise.venn(nrow(subset(df, BT.x == "Imbalanced Genes \n (according to replicate 1)")), 
#                      nrow(subset(df, BT.y == "Imbalanced Genes \n (according to replicate 1)")), 
#                      nrow(subset(df, BT.x == "Imbalanced Genes \n (according to replicate 1)" & BT.y == "Imbalanced Genes \n (according to replicate 1)")), 
#                      category = c("Tech Replicate #1", "Tech Replicate #2"), 
#                      cat.pos = c(-5, 5), cat.dist = rep(0.025, 2),
#                      lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
#                      scaled = TRUE)
# }
# # CC
# for(lib in unique(DF_forplot$BTBFCC$libprep)){
#   df = DF_forplot$BTBFCC[DF_forplot$BTBFCC$libprep == lib, ]
#   grid.newpage()
#   draw.pairwise.venn(nrow(subset(df, BT.x == "Imbalanced Genes \n (according to replicate 1)")), 
#                      nrow(subset(df, BT.y == "Imbalanced Genes \n (according to replicate 1)")), 
#                      nrow(subset(df, BT.x == "Imbalanced Genes \n (according to replicate 1)" & BT.y == "Imbalanced Genes \n (according to replicate 1)")), 
#                      category = c("Tech Replicate #1", "Tech Replicate #2"), 
#                      cat.pos = c(-5, 5), cat.dist = rep(0.025, 2),
#                      lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
#                      scaled = TRUE)
# }


library(venneuler)

plot(venneuler(na.omit(df[, c(2,4)])))
plot(venneuler(na.omit(df[, c(3,5)])))
plot(venneuler(na.omit(df[, c(2,4,6)])))
plot(venneuler(na.omit(df[, c(3,5,7)])))
plot(venneuler(na.omit(df[, c(2,4,6,8)])))
plot(venneuler(na.omit(df[, c(3,5,7,9)])))

plot(venneuler(na.omit(df[, c(2,4,6,8,10,12)])))
plot(venneuler(na.omit(df[, c(3,5,7,9,11,13)])))

par(mfrow=c(1, 3))

bt1 = venneuler(na.omit(df_NEB_SM10_SM100[, c(2,4,6,8,10,12)]))
bt2 = venneuler(na.omit(df_NEB_SM10_SM100[, 12+c(2,4,6,8,10,12)]))
bt3 = venneuler(na.omit(df_NEB_SM10_SM100[, 24+c(2,4,6,8,10,12)]))
bt1$labels <- rep("", 6)
bt2$labels <- rep("", 6)
bt3$labels <- rep("", 6)
plot(bt1)
plot(bt2)
plot(bt3)

btcc1 = venneuler(na.omit(df_NEB_SM10_SM100[, c(3,5,7,9,11,13)]))
btcc2 = venneuler(na.omit(df_NEB_SM10_SM100[, 12+c(3,5,7,9,11,13)]))
btcc3 = venneuler(na.omit(df_NEB_SM10_SM100[, 24+c(3,5,7,9,11,13)]))
btcc1$labels <- rep("", 6)
btcc2$labels <- rep("", 6)
btcc3$labels <- rep("", 6)
plot(btcc1)
plot(btcc2)
plot(btcc3)


par(mfrow=c(2, 3))

plot(bt1)
plot(bt2)
plot(bt3)
plot(btcc1)
plot(btcc2)
plot(btcc3)


par(mfrow=c(1,2))

bt_all = venneuler(na.omit(df_NEB_SM10_SM100[, c(c(2,4,6,8,10,12), 12+c(2,4,6,8,10,12), 24+c(2,4,6,8,10,12))]))
bt_all$labels <- rep("", 18)
plot(bt_all)
btcc_all = venneuler(na.omit(df_NEB_SM10_SM100[, c(c(3,5,7,9,11,13), 12+c(3,5,7,9,11,13), 24+c(3,5,7,9,11,13))]))
btcc_all$labels <- rep("", 18)
plot(btcc_all)

par(mfrow=c(1,2))

bt_all3 = venneuler(na.omit(df_NEB_SM10_SM100[, c(c(8,10,12), 12+c(8,10,12), 24+c(8,10,12))]))
bt_all3$labels <- rep("", 9)
plot(bt_all3)

btcc_all3 = venneuler(na.omit(df_NEB_SM10_SM100[, c(c(9,11,13), 12+c(9,11,13), 24+c(9,11,13))]))
btcc_all3$labels <- rep("", 9)
plot(btcc_all3)

par(mfrow=c(1,2))

bt_all2 = venneuler(na.omit(df_NEB_SM10_SM100[, c(c(4,10), 12+c(4,10), 24+c(4,10))]))
bt_all2$labels <- rep("", 6)
plot(bt_all2)

btcc_all2 = venneuler(na.omit(df_NEB_SM10_SM100[, c(c(5,11), 12+c(5,11), 24+c(5,11))]))
btcc_all2$labels <- rep("", 6)
plot(btcc_all2)

bt_all1 = venneuler(na.omit(df_NEB_SM10_SM100[, c(c(4), 12+c(4), 24+c(4))]))
bt_all1$labels <- rep("", 3)
plot(bt_all1)

btcc_all1 = venneuler(na.omit(df_NEB_SM10_SM100[, c(c(5), 12+c(5), 24+c(5))]))
btcc_all1$labels <- rep("", 3)
plot(btcc_all1)


par(mfrow=c(3,2))

bt1_2 = venneuler(na.omit(df_NEB_SM10_SM100[, c(4,10)]))
bt1_2$labels <- rep("", 2)
plot(bt1_2)
btcc1_2 = venneuler(na.omit(df_NEB_SM10_SM100[, c(5,11)]))
btcc1_2$labels <- rep("", 2)
plot(btcc1_2)

bt2_2 = venneuler(na.omit(df_NEB_SM10_SM100[, 12+c(4,10)]))
bt2_2$labels <- rep("", 2)
plot(bt2_2)
btcc2_2 = venneuler(na.omit(df_NEB_SM10_SM100[, 12+c(5,11)]))
btcc2_2$labels <- rep("", 2)
plot(btcc2_2)

bt3_2 = venneuler(na.omit(df_NEB_SM10_SM100[, 24+c(4,10)]))
bt3_2$labels <- rep("", 2)
plot(bt3_2)
btcc3_2 = venneuler(na.omit(df_NEB_SM10_SM100[, 24+c(5,11)]))
btcc3_2$labels <- rep("", 2)
plot(btcc3_2)


bt12 = venneuler(na.omit(df_NEB_SM10_SM100[, c(4,12+10)]))
bt12$labels <- rep("", 2)
plot(bt12)
btcc12 = venneuler(na.omit(df_NEB_SM10_SM100[, c(5,12+11)]))
btcc12$labels <- rep("", 2)
plot(btcc12)

bt13 = venneuler(na.omit(df_NEB_SM10_SM100[, c(4,24+10)]))
bt13$labels <- rep("", 2)
plot(bt13)
btcc13 = venneuler(na.omit(df_NEB_SM10_SM100[, c(5,24+11)]))
btcc13$labels <- rep("", 2)
plot(btcc13)

bt23 = venneuler(na.omit(df_NEB_SM10_SM100[, c(12+4,24+10)]))
bt23$labels <- rep("", 2)
plot(bt23)
btcc23 = venneuler(na.omit(df_NEB_SM10_SM100[, c(12+5,24+11)]))
btcc23$labels <- rep("", 2)
plot(btcc23)


par(mfrow=c(1,2))

bt_all5 = venneuler(na.omit(df_NEB_SM10_SM100[, c(c(4,6,8,10,12), 12+c(2,4,6,8,10,12), 24+c(2,4,6,8,10,12))]))
bt_all5$labels <- rep("", 18)
plot(bt_all5)
btcc_all5 = venneuler(na.omit(df_NEB_SM10_SM100[, c(c(5,7,9,11,13), 12+c(2,5,7,9,11,13), 24+c(2,5,7,9,11,13))]))
btcc_all5$labels <- rep("", 18)
plot(btcc_all5)


df_NEB_SM10_SM100_2 = rbind(df_NEB_SM10_SM100, df_NEB_SM10_SM100, df_NEB_SM10_SM100)
df_NEB_SM10_SM100_2[1:nrow(df_NEB_SM10_SM100), c(12+(2:13), 24+(2:13))] = FALSE
df_NEB_SM10_SM100_2[(nrow(df_NEB_SM10_SM100)+1):(2*nrow(df_NEB_SM10_SM100)), c(2:13, 24+(2:13))] = FALSE
df_NEB_SM10_SM100_2[(nrow(df_NEB_SM10_SM100)*2+1):(3*nrow(df_NEB_SM10_SM100)), c(2:13, 12+(2:13))] = FALSE


par(mfrow=c(1,1))

bt_btcc = venneuler(na.omit(df_NEB_SM10_SM100_2[, c(4:11, c(4:11)+12, c(4:11)+24)]))
bt_btcc$labels <- rep("", 18)
plot(bt_btcc)

# 
# HIST of (in 1),(in 2),(in 3),(in 4) ...
# 


