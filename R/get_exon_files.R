# About versions http://useast.ensembl.org/info/website/archives/assembly.html
library(biomaRt)
library("dplyr")
#mouse
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
exons_GRCm38 <- getBM(c("ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", "strand", "start_position", "end_position",
                        "transcript_start", "transcript_end",
                        "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end", 
                        "gene_biotype", "external_gene_name"), mart=ensembl)
exons_GRCm38 <- filter(exons_GRCm38, chromosome_name %in% c(1:19, "X", "Y"))
exons_GRCm38 <- filter(exons_GRCm38, gene_biotype %in% c("protein_coding", "lincRNA"))
exons_GRCm38$unique_id <- paste(exons_GRCm38$ensembl_transcript_id, exons_GRCm38$ensembl_exon_id, sep="_")
exons_GRCm38$chr <- paste0("chr", exons_GRCm38$chromosome_name)
exons_GRCm38 <- exons_GRCm38[order(exons_GRCm38$chromosome_name, exons_GRCm38$start_position, exons_GRCm38$ensembl_gene_id),]
write.table(exons_GRCm38, file="~/Dropbox (Partners HealthCare)/ASE/code/GRCm38.v96.exons_full.txt", quote = F, row.names = F, col.names = T, sep="\t")
exons_GRCm38_bed <- exons_GRCm38[,c('chr', 'exon_chrom_start', 'exon_chrom_end', 'unique_id')]
write.table(exons_GRCm38_bed, file="~/Dropbox (Partners HealthCare)/ASE/code/GRCm38.v96.exons.bed", quote = F, row.names = F, col.names = F, sep="\t")
#human
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
exons_GRCh38 <- getBM(c("ensembl_gene_id", "ensembl_transcript_id", "chromosome_name", "strand", "start_position", "end_position",
                        "transcript_start", "transcript_end",
                        "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end", 
                        "gene_biotype", "external_gene_name"), mart=ensembl)
exons_GRCh38 <- filter(exons_GRCh38, chromosome_name %in% c(1:19, "X", "Y"))
exons_GRCh38 <- filter(exons_GRCh38, gene_biotype %in% c("protein_coding", "lincRNA"))
exons_GRCh38$unique_id <- paste(exons_GRCh38$ensembl_transcript_id, exons_GRCh38$ensembl_exon_id, sep="_")
exons_GRCh38$chr <- paste0("chr", exons_GRCh38$chromosome_name)
exons_GRCh38 <- exons_GRCh38[order(exons_GRCh38$chromosome_name, exons_GRCh38$start_position, exons_GRCh38$ensembl_gene_id),]
write.table(exons_GRCh38, file="GRCh38.exons_full.txt", quote = F, row.names = F, col.names = T, sep="\t")
exons_GRCh38_bed <- exons_GRCh38[,c('chr', 'exon_chrom_start', 'exon_chrom_end', 'unique_id')]
write.table(exons_GRCh38_bed, file="GRCh38.exons.bed", quote = F, row.names = F, col.names = F, sep="\t")
