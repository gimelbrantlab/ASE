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



list_of_consts = list(mean(CC[[1]][6:15]), mean(CC[[2]]), mean(CC[[3]]))
list_of_libprepnames = list("NEBNext (100ng)", "SMARTseq (10ng)", "SMARTseq (0.1ng)")

df_NEB = CalculateTestWithHalf(data_30mln_noX[[1]], list_of_consts[[1]], list_of_libprepnames[[1]], thr=10)
df_SM10 = CalculateTestWithHalf(data_30mln_noX[[2]], list_of_consts[[2]], list_of_libprepnames[[2]], thr=10)
df_SM100 = CalculateTestWithHalf(data_30mln_noX[[3]], list_of_consts[[3]], list_of_libprepnames[[3]], thr=10)

df_NEB_SM10_SM100 = Reduce(function(x, y) merge(x, y, by="ID"), list(df_NEB, df_SM10, df_SM100))
dft = df_NEB_SM10_SM100
names(dft) = c("ID", paste0("rep", rep(1:18, each=2), '_', c("bt", "btCC")))
rownames(dft) = dft$ID


dft_2 = rbind(dft, dft, dft)
dft_2[1:nrow(dft), c(12+(2:13), 24+(2:13))] = FALSE
dft_2[(nrow(dft)+1):(2*nrow(dft)), c(2:13, 24+(2:13))] = FALSE
dft_2[(nrow(dft)*2+1):(3*nrow(dft)), c(2:13, 12+(2:13))] = FALSE
rownames(dft_2) = c(paste0(rownames(dft), ".1"), paste0(rownames(dft), ".2"), paste0(rownames(dft), ".3"))

eulij3exp = euler(na.omit(dft_2[, c('rep2_bt', 'rep2_btCC', 'rep4_bt', 'rep4_btCC',
                                    'rep8_bt', 'rep8_btCC', 'rep9_bt', 'rep9_btCC',
                                    'rep14_bt', 'rep14_btCC', 'rep15_bt', 'rep15_btCC'
)]),
shape='ellipse')
#eulij3exp_plt =
plot(eulij3exp,
     quantities = F,
     fills = F,
     edges = list(col=c("royalblue1","royalblue1","royalblue1","royalblue1",
                        "maroon2","maroon2","maroon2","maroon2",
                        "olivedrab3","olivedrab3","olivedrab3","olivedrab3"),
                  lty=c("dashed","solid","dashed","solid",
                        "dashed","solid","dashed","solid",
                        "dashed","solid","dashed","solid"),
                  lwd=3),
     labels = F,
     legend = F)
