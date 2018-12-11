exp16_6a <- read.delim("~/Dropbox (Partners HealthCare)/replicates_ASE/data/Sandberg_2014/GSM1278009_16cell_2pooled_split6a_expression.txt")
exp16_6b <- read.delim("~/Dropbox (Partners HealthCare)/replicates_ASE/data/Sandberg_2014/GSM1278010_16cell_2pooled_split6b_expression.txt")
exp16_8a <- read.delim("~/Dropbox (Partners HealthCare)/replicates_ASE/data/Sandberg_2014/GSM1278011_16cell_2pooled_split8a_expression.txt")
exp16_8b <- read.delim("~/Dropbox (Partners HealthCare)/replicates_ASE/data/Sandberg_2014/GSM1278012_16cell_2pooled_split8b_expression.txt")

exp16_6a <- exp16_6a[,c(1,5,6)]
colnames(exp16_6a) <- c("Gene_symbol", "rep1_ref", "rep1_alt")

exp16_6b <- exp16_6b[,c(1,5,6)]
colnames(exp16_6b) <- c("Gene_symbol", "rep2_ref", "rep2_alt")

geneCountTab <- merge(exp16_6a, exp16_6b, by.x="Gene_symbol", by.y="Gene_symbol")
colnames(geneCountTab) <- c("ensembl_gene_id","rep1_ref","rep1_alt","rep2_ref","rep2_alt")

dfDeltaAIPairwise = CreateMergedDeltaAIPairwiseDF(geneCountTab, what="exp_16")

DF = dfDeltaAIPairwise
DF$geneORsnp = "Gene"
DF$group = paste(DF$what, DF$geneORsnp, DF$ij)
head(DF)

DFQ = do.call(rbind, lapply(unique(DF$group),
                            function(r){
                              CreateObservedQuartilesDF(DF[DF$group == r, ],
                                                        P, ep=1.3, logbase=T,
                                                        coverageLimit=2000, group=r)
                            })
)
head(DFQ)

DFQ$method = sapply(as.character(DFQ$group), function(x){unlist(strsplit(x, '_'))[1]})
DFQ$ij = sapply(as.character(DFQ$group), function(x){paste(unlist(strsplit(x, ' '))[3:5], collapse = ' ')})
head(DFQ)

groupsIntercepts = lapply(unique(DFQ$Q), function(q){
  ddf = do.call(rbind, lapply(unique(DFQ$method), function(m){
    df = DFQ[DFQ$Q == q & DFQ$method == m, ]
    res = sapply(unique(df$ij), function(x){
      FitLmIntercept(df[df$ij == x, ], 40, morethan=10, logoutput=T)
    })
    names(res) = unique(df$ij)
    res
  }))
  row.names(ddf) = paste(unique(DFQ$method), "_ Q%:",
                         round(as.double(as.character(q)),4))
  2**ddf
})

do.call(rbind, lapply(groupsIntercepts, function(ddf){data.frame(round(2**ddf, 3))}))
