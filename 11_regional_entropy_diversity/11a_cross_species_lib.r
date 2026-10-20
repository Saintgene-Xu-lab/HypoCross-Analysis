
suppressMessages({
  # basic
  library(magrittr)
  library(Seurat)
  library(Matrix)
  library(future)
  library(furrr)
  library(patchwork)
  library(pals)
  library(future.apply)
  library(svglite)
  library(tidyverse)
  
  # specific lib
  library(RANN)
  library(igraph)
  library(harmony)
  library(RcppAnnoy)
  library(dtw)
  library(irlba)
})

get_nearby_embedding = function(seutb, k=200,
                                coord_labs=c("rx", "ry"), cluster_labs="seurat_clusters", smooth=T) {
    kk = k %/% 20 # number of neighbour to smooth
    seu = seutb
    seu[[cluster_labs]] = as.factor(seu[[cluster_labs]])
    spatial_coord = as.matrix(seu[, coord_labs])
    nn_res = nn2(data=spatial_coord, k=k, treetype="kd", searchtype="radius", radius=3000, eps=1)$nn.idx
    nn_res_smooth = nn_res[, 1:kk]


    empty_cell = as.character(nrow(seu)+1)
    mapt = c(seu[[cluster_labs]], factor(length(levels(seu[[cluster_labs]]))))
    names(mapt) = 1:(nrow(seu)+1)
    nres = replace(c(nn_res), c(nn_res)==0, nrow(seu)+1)
    res = mapt[nres]

    dim(res) = dim(nn_res)
    ress = table(rep(1:nrow(res), ncol(res)), c(res))
    # dim(res) = c(ncol(mseu), length(levels(mseu$seurat_clusters))+1)
    ress = ress[, 1:(ncol(ress)-1)]

    if (smooth) {
        ## smooth over `kk` neighbor
        ressm = apply(nn_res_smooth, 1, function(x) {
            colSums(ress[x, ])
        })
        ressm = t(ressm)
        ress = ressm
    }
    rownames(ress) = seu$cell
    ress
}

get_embedding_distance = function(embed, k=300, rsf=NULL) {
    set.seed(42)
    nccg = nn2(data=embed, k=k+1, treetype="kd", searchtype="priority", eps=1)
    weight = scale(nccg$nn.dists[,-1], center=F)
    weight = exp(-weight)

    i = rep(nccg$nn.idx[, 1], each=k)
    j = as.vector(t(nccg$nn.idx[, -1]))
    idx = j > i
    i = i[idx]
    j = j[idx]
    x = as.vector(x)[idx]

    wf = 1
    if (!is.null(rsf)) {
        wf = abs(log(-rsf[i]) - log(-rsf[j])) # weighted RSF factor
        wf = as.vector(10^(-wf))
        # wf = scale(wf)
        # sup = quantile(wf[wf>1], probs=0.99)
        # wf = ifelse(wf > sup, sup, wf)
    }
    x = x * wf
    
    
    # x = replace_na(x, 0)
    nccgsm = sparseMatrix(i=i,
                            j=j,
                            x=x,
                            symmetric=T)
    nccgsm = graph_from_adjacency_matrix(nccgsm, weighted=T, mode="undirected")
    nccgsm
}






















#’ calculate KNN embeding for cell 
#'
#' @param seutb Seurat object metadata coordinate
#' @param k number of neighbor
#' @param coord_labs vector of three element for slot in `seu`
#' @param cluster_labs labels for clustring
#' @return embedding
#' @noRd 
get_nearby_embedding_gs = function(seutb, k=100,
                                coord_labs=c("rx", "ry"), cluster_labs="seurat_clusters") {
    # kk = k %/% 20 # number of neighbour to smooth
    seu = seutb
    seu[[cluster_labs]] = as.factor(seu[[cluster_labs]])
    spatial_coord = as.matrix(seu[, coord_labs])
    nnr = nn2(data=spatial_coord, k=k, treetype="kd", searchtype="radius", radius=400, eps=1)
    nn_res = nnr$nn.idx
    # gsd = 1/mean(apply(t(nnr$nn.dists), 2, function(v) {
    #     sd(v[v!=0])
    # }))*100
    mask = nn_res != 0

    # nnr$nn.dists = 
    nd = nnr$nn.dists 
    nd[!mask] = 0 
    nd[mask] = dnorm(scale(nd[mask], center=F), sd=0.2)
    wm = apply(nd, 1, function(v) v/sum(v)*100)
    wm = t(wm)
    # nn_res_smooth = nn_res[, 1:kk]


    empty_cell = as.character(nrow(seu)+1)
    mapt = c(seu[[cluster_labs]], factor(length(levels(seu[[cluster_labs]]))))
    names(mapt) = 1:(nrow(seu)+1)
    nres = replace(c(nn_res), c(nn_res)==0, nrow(seu)+1)
    res = mapt[nres]
    dim(res) = dim(nn_res)

    ml = max(as.integer(levels(seu[[cluster_labs]])))
    hmat = t(sapply(1:nrow(res), function(i) {
        v=setNames(wm[i,], res[i,])
        v = xtabs(v~names(v))
        v[as.character(0:ml)]
    }))
    rownames(hmat) = seu$cell
    colnames(hmat) = as.character(0:ml)
    hmat[is.na(hmat)] = 0
    hmat = hmat[, 1:(ncol(ress)-1)]
    
    # ress = table(rep(1:nrow(res), ncol(res)), c(res))
    # dim(res) = c(ncol(mseu), length(levels(mseu$seurat_clusters))+1)
    # ress = ress[, 1:(ncol(ress)-1)]

    # if (smooth) {
    #     ## smooth over `kk` neighbor
    #     ressm = apply(nn_res_smooth, 1, function(x) {
    #         colSums(ress[x, ])
    #     })
    #     ressm = t(ressm)
    #     ress = ressm
    # }

    
    # rownames(ress) = seu$cell
    # ress
    hmat
}

guassian_kernel_seu = function(seu, k=15, coord=c("rx", "ry"), lambda=0.8) {
    set.seed(42)
    print(paste0("processing ", seu$chip[1]))
    metak = floor(k / 2)
    nnr = nn2(data=seu[[coord]], k=ifelse(ncol(seu) > k,k,ncol(seu)), 
            treetype="kd", searchtype="priority", 
            eps=1)
    idx = nnr$nn.idx
    dists = nnr$nn.dists
    wml = apply(dists, 1, function(v) {
        widx = v[metak]
        w = exp(-(v/widx)^2)
        w / sum(w)
    })
    seu = seu |> NormalizeData()
    mat = seu@assays$RNA@data
    agmat = sapply(1:nrow(idx), function(i) {
        as.vector(mat[, idx[i,]] %*% wml[, i])
        # rowSums(mat[,i[i!=0]])
    })
    colnames(agmat) = colnames(mat)
    rownames(agmat) = paste0("agg-", rownames(mat))
    # fmat = agmat
    fmat = rbind(scale(mat, center=F)*sqrt(1-lambda), scale(agmat, center=F)*sqrt(lambda))
    return(CreateSeuratObject(fmat, meta.data=seu@meta.data))
}
