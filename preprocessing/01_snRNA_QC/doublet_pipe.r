# /usr/bin/env R

### doublet finder pipe ###
 doublet_pipe = function(seu) {
  filter_rate = 0.075
  t = paramSweep_v3(seu, 1:10, sct=T, num.cores=1)
  tt = summarizeSweep(t, GT = FALSE)
  ttt = find.pK(tt)
  pk = as.numeric(levels(ttt[which(ttt$BCmetric == max(ttt$BCmetric)), "pK"]))[ttt[which(ttt$BCmetric == max(ttt$BCmetric)), "pK"]]
  pn = tt %>% filter(pK == pk)
  pn = as.numeric(levels(pn[which(pn$BCreal == max(pn$BCreal)), "pN"]))[pn[which(pn$BCreal == max(pn$BCreal)), "pN"]]
  mht = modelHomotypic(seu$seurat_clusters)
  nExp_poi = round(filter_rate*ncol(seu))
  nExp_poi.adj = round(nExp_poi*(1-mht))
  DFparam = paste(0.25, pk, nExp_poi, sep='_')
  nseu = doubletFinder_v3(seu, PCs=1:10, pN=0.25, pK=pk, nExp=nExp_poi, reuse.pANN=FALSE, sct=T)
  df_slot_name = paste("DF.classifications", DFparam, sep='_')
  df_slot_name2 = paste('pANN', DFparam, sep='_')
  nseu[["DFRes"]] = nseu[[df_slot_name]]
  nseu[["DFparam"]] = DFparam
  nseu[[df_slot_name]] = NULL
  nseu[[df_slot_name2]] = NULL
  nseu
}
