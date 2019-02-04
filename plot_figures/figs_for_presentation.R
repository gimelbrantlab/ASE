inDF18_copy <- inDF18[,c(1,14:25)]
colnames(inDF18_copy) <- colnames(inDF18)[1:13]

inDF18_copy$rep1_ai <- inDF18_copy$rep1_ref / (inDF18_copy$rep1_ref + inDF18_copy$rep1_alt)
inDF18_copy$rep2_ai <- inDF18_copy$rep2_ref / (inDF18_copy$rep2_ref + inDF18_copy$rep2_alt)
inDF18_copy$rep3_ai <- inDF18_copy$rep3_ref / (inDF18_copy$rep3_ref + inDF18_copy$rep3_alt)
inDF18_copy$rep4_ai <- inDF18_copy$rep4_ref / (inDF18_copy$rep4_ref + inDF18_copy$rep4_alt)
inDF18_copy$rep5_ai <- inDF18_copy$rep5_ref / (inDF18_copy$rep5_ref + inDF18_copy$rep5_alt)
inDF18_copy$rep6_ai <- inDF18_copy$rep6_ref / (inDF18_copy$rep6_ref + inDF18_copy$rep6_alt)

inDF18_copy$cov <- rowSums(inDF18_copy[, 2:13]) / 6

inDF18_copy <- inDF18_copy %>% select(ensembl_gene_id, rep1_ai, rep2_ai, rep3_ai,
                                      rep4_ai, rep5_ai, rep6_ai, cov)
inDF18_copy$sd <- apply(inDF18_copy[,2:7],1,sd)

ggplot(inDF18_copy, aes(x=cov, y=sd)) +
  geom_point() +
  theme_bw() +
  xlim(0, 30000)

v <- c(0,10,100,1000,10000)
inDF18_copy$cov_bin <- findInterval(inDF18_copy$cov, v)

g1 <- ggplot(inDF18_copy, aes(x=rep1_ai, y=rep2_ai)) +
  geom_point(size=0.5) +
  theme_bw() +
  coord_fixed() +
  xlab("AI: Technical Replicate 1") +
  ylab("AI: Technical Replicate 2")
save_plot(g1, file="corr_plot_bw.pdf", base_height = 4)

g2 <- ggplot(inDF18_copy, aes(x=rep1_ai, y=rep2_ai, col=6-cov_bin)) +
  geom_point(size=0.5) +
  theme_bw() +
  coord_fixed() +
  #theme(legend.position = "None") +
  xlab("AI: Technical Replicate 1") +
  ylab("AI: Technical Replicate 2")
save_plot(g2, file="corr_plot_col.pdf", base_height = 4)



df <- data.frame(group = rep(c("Above", "Below"), each=10), x = rep(1:10, 2), y = c(runif(10, 0, 1), runif(10, -1, 0)))
ggplot(df, aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity", position="identity")


SG2_1_genes <- read_delim("~/Dropbox (Partners HealthCare)/our_papers/5aza/data/SG2_1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
SG2_2_genes <- read_delim("~/Dropbox (Partners HealthCare)/our_papers/5aza/data/SG2_2.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
SG5_1_genes <- read_delim("~/Dropbox (Partners HealthCare)/our_papers/5aza/data/SG5_1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
SG5_2_genes <- read_delim("~/Dropbox (Partners HealthCare)/our_papers/5aza/data/SG5_2.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
SG2_1_genes <- SG2_1_genes %>% select(gene_id, TPM) %>% rename(TPM_SG2_1 = TPM)
SG2_2_genes <- SG2_2_genes %>% select(gene_id, TPM) %>% rename(TPM_SG2_2 = TPM)
SG5_1_genes <- SG5_1_genes %>% select(gene_id, TPM) %>% rename(TPM_SG5_1 = TPM)
SG5_2_genes <- SG5_2_genes %>% select(gene_id, TPM) %>% rename(TPM_SG5_2 = TPM)
TPMs <- full_join(SG2_1_genes, SG2_2_genes, "gene_id") %>% full_join(SG5_1_genes, "gene_id")  %>% full_join(SG5_2_genes, "gene_id")

df_compare_DE <- merge(df_compare, dd, by.x="ensembl_gene_id", by.y="gene_id")



TPMs_ai <- merge(TPMs, df_compare_DE, by.x="gene_id", by.y="ensembl_gene_id")
TPMs_ai$pval_sign <- ifelse(TPMs_ai$padj < 0.1, 1, 0)
TPMs_ai$col_gr <- ifelse(TPMs_ai$padj < 0.1, "DE", "NOTDE")
TPMs_ai$int <- ifelse(TPMs_ai$pval_sign == 0 & abs(TPMs_ai$ai_diff)>=0.25, 1, 0.6)
TPMs_ai$cov_av <- (TPMs_ai$exp2_ref + TPMs_ai$exp2_alt + TPMs_ai$exp5_ref + TPMs_ai$exp5_alt)/2

write.table(TPMs_ai, file="TPMs_ai_5_2.txt", quote = F, row.names = F)

group.colors <- c(DE = "red", NOTDE = "black")

ggplot(TPMs_ai, aes(x=exp2_ai, y=log2FoldChange, col=col_gr)) +
  geom_point(size=0.7) +
  theme_bw() +
  scale_color_manual(values = group.colors) +
  theme(legend.position = "None") +
  ylim(-2,2)



p1 <- ggplot(TPMs_ai, aes(x=ai_diff, y=log2FoldChange, col=col_gr)) +
  geom_point(size=0.7) +
  theme_bw() +
  theme(legend.position = "None") +
  scale_color_manual(values = group.colors) +
  ylim(-2,2) +
  xlim(-0.5, 0.5)
save_plot(p1, file="Fig_DE_vs_ai.pdf", base_aspect_ratio = 1.5)

p2 <- ggplot(TPMs_ai, aes(x=ai_diff, y=log2FoldChange, col=col_gr, alpha=int)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "None") +
  scale_color_manual(values = group.colors) +
  ylim(-2,2) +
  xlim(-0.5, 0.5)
save_plot(p2, file="Fig_DE_vs_ai_alpha.pdf", base_aspect_ratio = 1.5)

TPMs_ai_interest <- TPMs_ai[TPMs_ai$int == 1, c(15,10,11)]
TPMs_ai_interest$exp2_ai <- abs(TPMs_ai_interest$exp2_ai - 0.5)
TPMs_ai_interest$exp5_ai <- abs(TPMs_ai_interest$exp5_ai - 0.5)
TPMs_ai_interest$exp5_ai <- TPMs_ai_interest$exp5_ai * -1

TPMs_ai_interest_m <- melt(TPMs_ai_interest)
TPMs_ai_interest_m <- merge(TPMs_ai_interest_m, TPMs_ai[,c(18,15)], by.x="external_gene_name", by.y="external_gene_name")
v_bin <- c(10, 50, 100, 200, 500)
TPMs_ai_interest_m$cov_bin <- 0.25 * findInterval(TPMs_ai_interest_m$cov_av, v_bin) + 0.1

g1 <- ggplot(TPMs_ai_interest_m, aes(y=value, x=external_gene_name, fill=variable)) +
  geom_bar(stat="identity", position="identity") +
  theme(legend.position = "None", axis.text.x = element_text(angle = 90, vjust=0.5)) +
  ylab("Ratio") +
  ylim(-1,1)
  #+ coord_flip()
save_plot(g1, file="Fig_5_2_notDE.pdf", base_aspect_ratio = 2)

g2 <- ggplot(TPMs_ai_interest_m, aes(y=value, x=external_gene_name, fill=variable, alpha=cov_bin)) +
  geom_bar(stat="identity", position="identity") +
  theme(legend.position = "None", axis.text.x = element_text(angle = 90, vjust=0.5)) +
  ylab("Ratio") +
  ylim(-1,1) +
  geom_hline(yintercept = 0.5, linetype="dashed", col="gray") +
  geom_hline(yintercept = -0.5, linetype="dashed", col="gray")
save_plot(g2, file="Fig_5_2_notDE_alpha.pdf", base_aspect_ratio = 2)

g3 <- ggplot(TPMs_ai_interest_m, aes(y=value, x=external_gene_name, fill=variable, alpha=cov_bin)) +
  geom_bar(stat="identity", position="identity") +
  theme(legend.position = "None", axis.text.x = element_text(angle = 90, vjust=0.5)) +
  ylab("Ratio") +
  ylim(-0.5,0.5) +
  geom_hline(yintercept = 0.35, linetype="dashed", col="gray") +
  geom_hline(yintercept = -0.35, linetype="dashed", col="gray")
save_plot(g3, file="Fig_5_2_notDE_alpha_abs.pdf", base_aspect_ratio = 2)

## plot with genes in 5aza
t <- read.delim("genes_notDE_5aza.txt")
t$type <- as.factor(t$type)
levels(t$type) <- c("B", "A")
g5 <- ggplot(t, aes(y=AI, fill=type, x=num)) + geom_bar(stat="identity", position="dodge") + geom_errorbar(aes(ymin=AI-CI, ymax=AI+CI), width=.2,position=position_dodge(.9)) + coord_flip() + scale_x_discrete(limits = rev(levels(as.factor((t$num))))) + geom_hline(yintercept = 1)
save_plot(g5, file="Fig_genes_notDE_5aza.pdf")
