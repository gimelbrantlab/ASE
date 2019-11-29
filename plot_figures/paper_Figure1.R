setwd("~/Dropbox (Partners HealthCare)/replicates_ASE/code/ASE/plot_figures")


library(tidyverse)
library(fitdistrplus)
library(cowplot)
library(ggrepel)

source("../R/ASE_functions.R")
source("../R/PerformAIAnalysis_CC.v0.3.R")
source("../R/DownstreamASE.R")
source("../plot_figures/utilPlots.R")
source("../plot_figures/takedata30mln.R")
source("../plot_figures/paper_Figure1_c.R")


## Gendrel = GetGatkPipelineTabs("../../../data/Gendrel/Gendrel_52604442.1_processed_gene_extended2.txt", c(4), multiple = T)
## chrgenes <- read.delim("../../../data/kidney/submln/MLN30_SG3_N955_NEB_R1_merged_v2.mln30_trial5_processed_gene_extended2.txt")[, c("ensembl_gene_id", "chr")]
## Gendrel <- Gendrel[Gendrel$ensembl_gene_id %in% chrgenes[chrgenes$chr!="chrX" & chrgenes$chr!="chrY", "ensembl_gene_id"], ]


# Select pair of replicates for analysis (for seed=1 it is replicates 2 and 6)
# set.seed(1)
# sample2reps30mln = sample(0:5, 2, replace=F)*5 + sample(1:5, 2, replace=T)
# sample2reps30mln_6 = sort((sample2reps30mln-1) %/% 5 + 1)
# Fixed 2 and 3 rep, first samples
sample2reps30mln = c(6, 11)
sample2reps30mln_6 = c(2, 3)

set.seed(1)
CCn_sample2reps30mln = sum((5:1)[0:(sample2reps30mln_6[1]-1)]) + sample2reps30mln_6[2]-sample2reps30mln_6[1]

CC_sample2_30mln = lapply(CC, function(x){x[CCn_sample2reps30mln]})

data_30mln_noX_sample2 = lapply(data_XXmln_noX, function(x){
  x[, sort(c(1, sample2reps30mln*2, sample2reps30mln*2+1))]
})

## CC_Gendrel <- c(mean(1.554912,1.580067), mean(1.581116,1.547620))

# Get data frame with AIs and coverages for 2 selected replicates + tests outputs (T/F)

# Get data for NEB

DF_NEB_BT <- GetDataForExperiment_BT(df = data_30mln_noX_sample2, CC_df = CC_sample2_30mln, exp_name = "NEB", exp_n = 1)

## DF_Gend_BT <- GetDataForExperiment_BT(df = Gendrel, CC_df = CC_Gendrel, exp_name = "Gendrel", exp_n = 2)


# Plot figures with classification (middle panel)

figure_1D = ggplot(DF_NEB_BT[DF_NEB_BT$meanCOV.x < 9999, ], aes(meanCOV.x, AI.x)) +
  geom_point(aes(color=BT.x), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  facet_grid(~ BT.x) +
  theme_bw() +
  #ggtitle("Replicate 1") +
  labs(x = "Total Allelic Counts", y = "Gene AI - Replicate 1", color = "Binomial Test") +
  scale_x_continuous(trans='log10') +
  theme(legend.position="None",  text = element_text(size=18),
        strip.background = element_rect(fill="#E0E0E0"))

figure_1E = ggplot(DF_NEB_BT[DF_NEB_BT$meanCOV.y < 9999, ], aes(meanCOV.y, AI.y)) +
  geom_point(aes(color=BT.y), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  facet_grid(~ BT.x) +
  theme_bw() +
  #ggtitle("Replicate 2 relative to replicate 1") +
  labs(x = "Total Allelic Counts", y = "Gene AI - Replicate 2", color = "Binomial Test") +
  scale_x_continuous(trans='log10') +
  theme(legend.position="None", text = element_text(size=18),
        strip.background = element_rect(fill="#E0E0E0"))


# Repeat analysis for SMART-seq (right panel)

DF_SMART_10_BT <- GetDataForExperiment_BT(df = data_30mln_noX_sample2, CC_df = CC_sample2_30mln, exp_name = "SMARTseq 10ng", exp_n = 2)

figure_1F <- ggplot(DF_SMART_10_BT[DF_SMART_10_BT$meanCOV.y < 9999 & DF_SMART_10_BT$BT.x=="Imbalanced Genes \n (according to replicate 1)", ], aes(meanCOV.y, AI.y)) +
  geom_point(aes(color=BT.y), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  facet_grid(~ BT.x) +
  theme_bw() +
  labs(x = "Total Allelic Counts", y = "Gene AI  - Replicate 2", color = "Binomial Test") +
  scale_x_continuous(trans='log10') +
  theme(legend.position="None", text = element_text(size=18),
        strip.background = element_rect(fill="#E0E0E0"))

DF_SMART_0.1_BT <- GetDataForExperiment_BT(df = data_30mln_noX_sample2, CC_df = CC_sample2_30mln, exp_name = "SMARTseq 0.1ng", exp_n = 3)

figure_1G <- ggplot(DF_SMART_0.1_BT[DF_SMART_0.1_BT$meanCOV.y < 9999 & DF_SMART_0.1_BT$BT.x=="Imbalanced Genes \n (according to replicate 1)", ], aes(meanCOV.y, AI.y)) +
  geom_point(aes(color=BT.y), size=0.8) +
  scale_color_manual(values=c("salmon","royalblue1")) +
  facet_grid(~ BT.x) +
  theme_bw() +
  labs(x = "Total Allelic Counts", y = "Gene AI  - Replicate 2", color = "Binomial Test") +
  scale_x_continuous(trans='log10') +
  theme(legend.position="None", text = element_text(size=18),
        strip.background = element_rect(fill="#E0E0E0"))

# Left panel

DF_NEB_BT$COV = rowMeans(DF_NEB_BT[, c("meanCOV.x", "meanCOV.y")])

figure_1A = ggplot(DF_NEB_BT, aes(AI.x, AI.y)) +
  geom_point(aes(color=COV), size=0.8) +
  scale_color_gradient(low="lightgray", high="black", trans = "log2",
                       breaks = c(10, 500, 5000),
                       labels = c(10, 500, 5000)) +
  theme_bw() +
  coord_fixed() +
  labs(x = "Gene AI - Replicate 1", y = "Gene AI - Replicate 2", color = "Allelic Counts") +
  theme(legend.key = element_blank(),legend.background=element_blank(),
        legend.position=c(0.2,0.8), text = element_text(size=18))

figure_1B = ggplot(DF_NEB_BT, aes(meanCOV.x, meanCOV.y)) +
  geom_point(color="black", size=0.8) +
  theme_bw() +
  coord_fixed() +
  labs(x = "Total Allelic Counts - Replicate 1", y = "Total Allelic Counts - Replicate 2") +
  theme(legend.position="bottom", text = element_text(size=18)) +
  xlim(0,10000) + ylim(0,10000)


figure_1C_pre = plot_grid(plt_sim, NULL, plt_same, NULL, plt_diff, nrow=1, rel_widths = c(1,0.3,1.05,0.3,1.2))
figure_1C = plot_grid(NULL, ggdraw(figure_1C_pre) +
                        draw_label(NumP[1,1], size = 18, fontface = 'bold', x = (1/2)/(sum(c(1,0.3,1.05,0.3,1.2))), y = 0.53) +
                        draw_label(NumP[1,2], size = 18, fontface = 'bold', x = (1/2)/(sum(c(1,0.3,1.05,0.3,1.2))), y = 0.47) +
                        draw_label(NumP[2,1], size = 18, fontface = 'bold', x = (1+0.3+1.05/2)/(sum(c(1,0.3,1.05,0.3,1.2))), y = 0.53) +
                        draw_label(NumP[2,2], size = 18, fontface = 'bold', x = (1+0.3+1.05/2)/(sum(c(1,0.3,1.05,0.3,1.2))), y = 0.47) +
                        draw_label(NumP[3,1], size = 18, fontface = 'bold', x = (1+0.3+1.05+0.3+1.2/2)/(sum(c(1,0.3,1.05,0.3,1.2))), y = 0.53) +
                        draw_label(NumP[3,2], size = 18, fontface = 'bold', x = (1+0.3+1.05+0.3+1.2/2)/(sum(c(1,0.3,1.05,0.3,1.2))), y = 0.47),
                      ncol = 1, rel_heights = c(1,2)
)

PLT_fig1 = plot_grid(
  plot_grid(plot_grid(figure_1A, figure_1B,
                      rel_widths = c(1,1), align = 'v', scale = c(0.9, 0.9), nrow = 1, labels = c("A", "B")),
            plot_grid(figure_1C, scale = c(0.9), nrow = 1, labels = c("C")),
            ncol = 1, rel_widths = c(1,1)),
  plot_grid(figure_1D, figure_1F, figure_1E, figure_1G,
            rel_widths = c(1,0.55,1,0.55), align = 'h', scale = c(0.9, 0.9, 0.9, 0.9), nrow = 2, labels = c("D", "F", "E", "G")),
  align = 'vh', rel_widths = c(2, 1.55)
)
# PLT_fig1 = plot_grid(
#   plot_grid(figure_1A, figure_1B, figure_1D, figure_1F,
#             rel_widths = c(1,1,1,0.55), align = 'h', scale = c(0.9, 0.9, 0.9, 0.9), nrow = 1, labels = c("A", "B", "D", "F")),
#   plot_grid(figure_1C, figure_1E, figure_1G,
#             rel_widths = c(2,1,0.55), align = 'v', scale = c(0.9, 0.9, 0.9, 0.9), nrow = 1, labels = c("C","E", "G")),
#   nrow = 2, align = 'vh'
# )
PLT_fig1

cowplot::save_plot(PLT_fig1, file="fig.1_v2.pdf", base_height = 12, base_width = 24)



