suppressMessages(library(stringr))

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
    pos=read.csv(label) %>%
            dplyr::rename(area_name=Acronym)
    # pos = pos[str_detect(pos$origin_name, "HyTh"),] %>%
	# 	dplyr::select(global_id, name) %>%
	# 	rename(area_name=name, gene_area=global_id) %>% 
	# 	dplyr::distinct()
    # pos$area_name = sapply(str_split(pos$area_name, '-'), function(x) {x[1]})
  
    # read.delim("/home/majy01/data/cell2loc/HYRegions_All.tsv")
        
    spatial_file=paste0(rawdbDir,'/',f)
    spatial=read.csv(file=spatial_file)

    # for HYPO data, replace gene_area with RegionID 
    spatial %<>% mutate(cell_id=cell_label)
    allname=levels(as.factor(pos$area_name))

    stemname=allname
    pos=pos[pos$area_name %in% stemname,]
    spatial=left_join(spatial,pos)

    # skip NA filter for data with no regionID
    spatial=na.omit(spatial)
    # spatial$area_name = ifelse(is.na(spatial$area_name), "unlabel", spatial$area_name)
    if(nrow(spatial)!=0){
    secRawTb=spatial
    secDnbSumTb = secRawTb %>% group_by(x, y) %>% summarize(
                UMICount = sum(umi_count),
                detectGeneCount = length(unique(gene)),
                label = cell_id[1],
                areaID = area_name[1],
                rx = mean(rx),
                ry = mean(ry)
            ) %>% ungroup() # Be careful to ungroup before mutate
            # Add a DNB density information
            secDnbSumTb = secDnbSumTb %>% mutate(
                dnbScaledDensity = getScaledDensity(x, y, n=100)
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
            x = as.integer(mean(x)), # Use DNB mean position as cell position
            y = as.integer(mean(y)),
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
    secCellSumTb = secCellTb %>% group_by(cell_id, secCellId, x, y, rx, ry) %>% summarize(
            detectGeneCount = n(),
            UMICount = sum(UMICount	),
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
        dimName = list(levels(secCellTb$gene), levels(secCellTb$secCellId))
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
    cellBinSeurat[["mt"]] = 0 #PercentageFeatureSet(cellBinSeurat, pattern="^mt-")
    cellBinSeurat[["cf_ratio"]] = cellBinSeurat$nCount_RNA / cellBinSeurat$nFeature_RNA
    cname = str_extract(f, "T[[:number:]]+")
    cellBinSeurat$chip = cname
    filname=gsub('txt','rds',f)
    signame=f
    outFile = paste0(outCountMxDir, "/", cname, ".rds")
    outSig = paste0(outCountMxDir, "/", signame)

    saveRDS(cellBinSeurat, file=outFile)
    file.create(outSig)
    }
}
args = commandArgs(trailingOnly=T)
rawdbDir=args[1]
outCountMxDir=args[2]
label=args[3]

rawdbFile = basename(Sys.glob(paste0(rawdbDir, "/*.csv")))
unDoneFile = rawdbFile
DoneFile = basename(Sys.glob(paste0(outCountMxDir, "/*.csv")))
unDoneFile = unDoneFile[unDoneFile %in% setdiff(unDoneFile, DoneFile)]

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

# reverse chip order accroding to AP order
mseu$ori_chip = str_split_i(mseu$chip, '_', 2)

saveRDS(mseu, paste0(outCountMxDir, "/mseu.rds"), compress="gzip")

