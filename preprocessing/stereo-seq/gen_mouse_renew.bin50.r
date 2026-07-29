
suppressMessages(library(stringr))
# suppressMessages(library(viridis))
# suppressMessages(library(cowplot))
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(parallel)
library(magrittr)
library(MASS)
library(stringr)
library(Matrix)
library(future)
library(future.apply)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

getScaledDensity = function(x, y, ...) {
    if (length(x) < 2 | length(y) < 2) {
        relDens = rep(0, length(x))
    } else {
        dens = get_density(x, y, ...)
        relDens = as.vector(scale(dens))
    }
    return(relDens)
}



generateBinCount = function(f, rawdbDir, outCountMxDir, label){
    message(f)
    pos=read.csv(label) %>% rename(area_name=Acronym)
        
    spatial_file=paste0(rawdbDir,'/',f)
    spatial=read.csv(file=spatial_file) |> mutate(cell_id = paste(x, y, sep='_'))

    # for HYPO data, replace gene_area with RegionID 
    #spatial %<>% mutate(gene_area=RegionID)
    allname=levels(as.factor(pos$area_name))

    # stemname=allname
    # pos=pos[pos$area_name %in% stemname,]
    spatial=left_join(spatial,pos)
    #spatial=na.omit(spatial)
    spatial=spatial[!is.na(spatial$area_name),]

    if(nrow(spatial)!=0){
    secRawTb=spatial
    secDnbSumTb = secRawTb %>% group_by(rx, ry) %>% summarize(
        UMICount = sum(umi_count),
        detectGeneCount = length(unique(gene)),
        label = cell_id[1],
        areaID = area_name[1],
        rx = mean(rx),
        ry = mean(ry)
    ) %>% ungroup() # Be careful to ungroup before mutate
            # Add a DNB density information
    secDnbSumTb = secDnbSumTb %>% mutate(
        dnbScaledDensity = getScaledDensity(rx, ry, n=100)
    )
            #write_tsv(secDnbSumTb, outDnbSumFile)
            # all(secDnbSumTb$cellNuniqueCount) == 1 # Should be TRUE
            # all(secDnbSumTb$areaIDuniqueCount) == 1 # Should be TRUE
    # Plot check files
    secDnbSumTb$labelFlag = ifelse(secDnbSumTb$label != 0, "cellDetected", "nonCell")
    set.seed(123)
# For each cell, we use the areaID with most DNB support.
    secCellSuppTb = secDnbSumTb %>% group_by(label) %>% summarize(
        mainAreaID = ifelse(
            all(is.na(areaID)), NA,
            names(sort(table(areaID), decreasing=T))[1]
        ),
        rx = as.integer(mean(rx)),
        ry = as.integer(mean(ry))
    ) %>% ungroup()
# Summary each cell gene count.
    secCellTb = secRawTb %>% group_by(gene, cell_id) %>% summarize(
        UMICount = sum(umi_count)
    ) %>% ungroup()
    colnames(secCellSuppTb)[1]='cell_id'
    # Remove cell_id 0
    secCellTb = subset(secCellTb, cell_id != 0)
    secCellTb = secCellTb %>% left_join(secCellSuppTb)
    secCellTb = secCellTb %>% mutate(
        secCellId = paste(cell_id)
    )
    secCellSumTb = secCellTb %>% group_by(cell_id, secCellId, rx, ry) %>% summarize(
        detectGeneCount = n(),
        UMICount = sum(UMICount),
        mainAreaID = mainAreaID[1]
    ) %>% ungroup()
        #outCellSumFile = paste0(outCellSumDir, "/", sub("txt", "spatialCellSum.txt", f))
        #write_tsv(secCellSumTb, outCellSumFile)
        # Plot check files
        #png(outFigFile, height=1000, width=1000, type="cairo")
        #print(gp1)
        #dev.off()
    secCellTb = secCellTb %>% mutate(
        geneID = as.factor(gene),
        secCellId = as.factor(secCellId)
    )
    # dimName = list(levels(secCellTb$gene), levels(secCellTb$secCellId))
    secCellTb$gene=as.factor(secCellTb$gene)  
    dimName = list(levels(secCellTb$gene), levels(secCellTb$secCellId))
    secCellWideCountMx = sparseMatrix(
        as.integer(secCellTb$gene), as.integer(secCellTb$secCellId), 
        x=secCellTb$UMICount, dimnames=dimName
    )
        #outWideCountMxFile = paste0(outWideCountMxDir, "/", sub("txt", "spatialCell.wideCountMx.rda", f))
        #save(secCellWideCountMx, file=outWideCountMxFile)
        #print(paste("Done spatial cell sum on", secRawFile))
    cellBinSeurat = CreateSeuratObject(counts=secCellWideCountMx)
    secCellSumTb$cell_id=as.character(secCellSumTb$cell_id)
    secCellSumTb=secCellSumTb[order(secCellSumTb$cell_id),]
    cellBinSeurat@meta.data$x=secCellSumTb$x
    cellBinSeurat@meta.data$y=secCellSumTb$y
    cellBinSeurat@meta.data$rx=secCellSumTb$rx
    cellBinSeurat@meta.data$ry=secCellSumTb$ry
    cellBinSeurat@meta.data$mainAreaID=secCellSumTb$mainAreaID
    cellBinSeurat[["mt"]] = PercentageFeatureSet(cellBinSeurat, pattern="^mt-")
    cellBinSeurat[["cf_ratio"]] = cellBinSeurat$nCount_RNA / cellBinSeurat$nFeature_RNA
    #cellBinSeurat_qc=subset(cellBinSeurat,subset = nFeature_RNA > 100 & nCount_RNA > 1.5 & percent.mt<15)
    #cellBinSeurat_qc
    #saveRDS(cellBinSeurat_qc,"/OceanStor100D/home/sunyidi_lab/liuyuxuan/spa_seurat/t317/t317_all.rds")
    cname = str_extract(f, "T[[:number:]]+")
    cellBinSeurat$chip = cname
    filname=gsub('csv','rds',f)
    signame=f
    outFile = paste0(outCountMxDir, "/", cname, ".rds")
    outSig = paste0(outCountMxDir, "/", signame)

    saveRDS(cellBinSeurat, file=outFile)
    file.create(outSig)
    #saveRDS(cellBinSeurat_qc,"/OceanStor100D/home/sunyidi_lab/liuyuxuan/spa_seurat/t317/piceng_fil.rds")
    }
}
args = commandArgs(trailingOnly=T)
#f='total_gene_T291_mouse_f001_2D_mouse1-20220823.txt'
rawdbDir=args[1]
#rawdbDir='/home/liuyx/liuyuxuan/spa_seurat/823_new_spa/all_mop_rds'
outCountMxDir=args[2]
label=args[3]

rawdbFile = basename(Sys.glob(paste0(rawdbDir, "/*.csv")))
unDoneFile = rawdbFile
DoneFile = basename(Sys.glob(paste0(outCountMxDir, "/*.csv")))
unDoneFile = unDoneFile[unDoneFile %in% setdiff(unDoneFile, DoneFile)]
# unDoneFile = unDoneFile[unDoneFile %in% setdiff(unDoneFile[49:83], DoneFile)]

message("unDone: ", unDoneFile)

thread = 4
options(future.globals.maxSize = +Inf)
plan(multicore, workers=thread)
future_lapply(
    unDoneFile, 
    function(f) {generateBinCount(f, rawdbDir, outCountMxDir, label)}
)
rdsl = Sys.glob(paste0(outCountMxDir, "/*.rds"))
rdsl = lapply(rdsl, readRDS)
mseu = merge(rdsl[[1]], rdsl[-1])
saveRDS(mseu, paste0(outCountMxDir, "/mseu.rds"), compress="gzip")
