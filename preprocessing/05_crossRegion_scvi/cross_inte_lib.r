
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






# ic <- function(harmonyObj, cluster_prior=NULL) {
#     if (harmonyObj$ran_setup == FALSE) {
#         stop('before initializing cluster, run setup')
#     }
#     if (!is.null(cluster_prior)) {
#         if (ncol(cluster_prior) != harmonyObj$N) {
#             stop('cluster_prior must be defined by N cells')
#         }
#         if (nrow(cluster_prior) > harmonyObj$K) {
#             stop('cluster_prior cannot contain more than K clusters')
#         }
#         C <- nrow(cluster_prior)
#         harmonyObj$Y <- matrix(0, harmonyObj$d, harmonyObj$K)
#         harmonyObj$Y[, seq_len(C)] <- compute_Y(harmonyObj$Z_cos,
#                                                            cluster_prior)
#         harmonyObj$R <- matrix(0, harmonyObj$K, harmonyObj$N)
#         harmonyObj$R[seq_len(nrow(cluster_prior)), ] <- cluster_prior
        

#         ## if needed, initialize K-C clusters        
#         if (C < harmonyObj$K) {
#             Ynew <- t(stats::kmeans(t(harmonyObj$Z_cos), 
#                                     centers = harmonyObj$K - C,
#                                     iter.max = 250, nstart = 100)$centers)
#             harmonyObj$Y[, seq(1+C, harmonyObj$K)] <- Ynew
#         }
        
#         harmonyObj$init_cluster_cpp(C)
#     } else {
#         harmonyObj$Y <- t(stats::kmeans(t(harmonyObj$Z_cos), 
#                                         centers = harmonyObj$K,
#                                     iter.max = 250, nstart = 100)$centers)
#         harmonyObj$init_cluster_cpp(0)
#     }

# }
# environment(ic) <- environment(harmony:::init_cluster)


hm = function (data_mat, meta_data, vars_use, do_pca = TRUE, npcs = 20, 
    theta = NULL, lambda = NULL, sigma = 0.1, nclust = NULL, 
    tau = 0, block.size = 0.05, max.iter.harmony = 10, max.iter.cluster = 200, 
    epsilon.cluster = 1e-05, epsilon.harmony = 1e-04, plot_convergence = FALSE, 
    return_object = FALSE, verbose = TRUE, reference_values = NULL, 
    cluster_prior=NULL, ...) 
{
    if (!(is(meta_data, "data.frame") | is(meta_data, "DataFrame"))) {
        if (length(meta_data) %in% dim(data_mat)) {
            meta_data <- data.frame(batch_variable = meta_data)
            vars_use <- "batch_variable"
        }
        else {
            stop("meta_data must be either a data.frame or a vector with batch \n                values for each cell")
        }
    }
    if (is.null(vars_use) | any(!vars_use %in% colnames(meta_data))) {
        msg <- gettextf("must provide variables names (e.g. vars_use=%s)", 
            sQuote("stim"))
        stop(msg)
    }
    if (do_pca) {
        if (ncol(data_mat) != nrow(meta_data)) {
            data_mat <- Matrix::t(data_mat)
        }
        pca_res <- data_mat %>% harmony:::scaleData() %>% irlba::prcomp_irlba(n = npcs,  maxit = 1000,
        # work = nv + 20,
            retx = TRUE, center = FALSE, scale. = FALSE)
        data_mat <- pca_res$rotation %*% diag(pca_res$sdev)
    }
    N <- nrow(meta_data)
    cells_as_cols <- TRUE
    if (ncol(data_mat) != N) {
        if (nrow(data_mat) == N) {
            data_mat <- t(data_mat)
            cells_as_cols <- FALSE
        }
        else {
            stop("number of labels do not correspond to number of \n                samples in data matrix")
        }
    }
    if (is.null(nclust)) {
        nclust <- min(round(N/30), 100)
    }
    if (is.null(theta)) {
        theta <- rep(2, length(vars_use))
    }
    else if (length(theta) != length(vars_use)) {
        stop("Please specify theta for each variable")
    }
    if (is.null(lambda)) {
        lambda <- rep(1, length(vars_use))
    }
    else if (length(lambda) != length(vars_use)) {
        stop("Please specify lambda for each variable")
    }
    if (length(sigma) == 1 & nclust > 1) {
        sigma <- rep(sigma, nclust)
    }
    phi <- Reduce(rbind, lapply(vars_use, function(var_use) {
        t(harmony:::onehot(meta_data[[var_use]]))
    }))
    N_b <- rowSums(phi)
    Pr_b <- N_b/N
    B_vec <- Reduce(c, lapply(vars_use, function(var_use) {
        length(unique(meta_data[[var_use]]))
    }))
    theta <- Reduce(c, lapply(seq_len(length(B_vec)), function(b) rep(theta[b], 
        B_vec[b])))
    theta <- theta * (1 - exp(-(N_b/(nclust * tau))^2))
    lambda <- Reduce(c, lapply(seq_len(length(B_vec)), function(b) rep(lambda[b], 
        B_vec[b])))
    lambda_mat <- diag(c(0, lambda))
    if (!is.null(reference_values)) {
        idx <- which(row.names(phi) %in% reference_values)
        cells_ref <- which(colSums(phi[idx, , drop = FALSE] == 
            1) >= 1)
        b_keep <- which(!row.names(phi) %in% reference_values)
        phi_moe <- phi[b_keep, , drop = FALSE]
        phi_moe[, cells_ref] <- 0
        phi_moe <- rbind(rep(1, N), phi_moe)
        lambda_mat <- lambda_mat[c(1, b_keep + 1), c(1, b_keep + 
            1)]
    }
    else {
        phi_moe <- rbind(rep(1, N), phi)
    }
    harmonyObj <- new(harmony, 0)
    harmonyObj$setup(data_mat, phi, phi_moe, Pr_b, sigma, theta, 
        max.iter.cluster, epsilon.cluster, epsilon.harmony, nclust, 
        tau, block.size, lambda_mat, verbose)
    init_cluster(harmonyObj, cluster_prior) #ic
    harmonize(harmonyObj, max.iter.harmony, verbose)
    if (plot_convergence) 
        graphics::plot(HarmonyConvergencePlot(harmonyObj))
    if (return_object) {
        return(harmonyObj)
    }
    else {
        res <- as.matrix(harmonyObj$Z_corr)
        row.names(res) <- row.names(data_mat)
        colnames(res) <- colnames(data_mat)
        if (!cells_as_cols) 
            res <- t(res)
        return(res)
    }
}

environment(hm) <- environment(HarmonyMatrix)



#’ calculate KNN embeding distance with harmony 
#'
#' @param embed KNN embedding from `get_nearby_embedding`
#' @param k number of neighbor
#' @param npcs pca reduaction of original dataset
#' @param get_ebd get embedding (default false)
#' @param batch labels for batch correction
#' @return graph object for clustering
#' @noRd 
get_harmony_graph = function(embed, k=500, npcs=30, get_ebd=F, batch) {
    hres = matrix(c(embed), ncol=ncol(embed))
    # hres = t(apply(hres, 1, scale))
    rownames(hres) = rownames(embed)
    set.seed(42)
    # approxi_pca = irlba(scale(hres), nv=npcs)
    # sdev <- approxi_pca$d/sqrt(max(1, ncol(hres) - 1))

    hebd <- HarmonyMatrix(embed,#approxi_pca$u %*% diag(approxi_pca$d), 
                            tibble(batch=batch, rn=rownames(hres)) %>% column_to_rownames("rn"), 
                            vars_use="batch", do_pca=T,
                            verbose=T)

    if (get_ebd) return(hebd)

    hebd_annoy <- new(AnnoyEuclidean, ncol(hebd))
    hebd_annoy$setSeed(42)
    hebd_annoy$setVerbose(1)
    for (i in 1:nrow(hebd)) {
        hebd_annoy$addItem(i-1, hebd[i,])
    }
    hebd_annoy$build(100)
    nn = future_sapply(1:nrow(hebd), function(i) {
        hebd_annoy$getNNsByVectorList(hebd[i,], k, -1, T)
    })
    hebd_annoy$unload()
    x = rep(1:nrow(hebd), each=k)
    y = unlist(sapply(1:nrow(hebd), function(i) {nn[[2*i-1]]}, simplify=F))
    v = unlist(sapply(1:nrow(hebd), function(i) {nn[[2*i]]}, simplify=F))
    # v = lapply(v, function(vv) {exp(-vv)})
    v = exp(-v)
    idx = (y+1) < x
    nccgsm = sparseMatrix(i=x[idx],
                        j=y[idx] + 1,
                        x=v[idx],
                        symmetric=T)
    nccgsm = graph_from_adjacency_matrix(nccgsm, weighted=T, mode="undirected")
    nccgsm
    # leid = cluster_leiden(nccgsm, objective_function="CPM",
    #                         resolution_parameter=0.01, n_iterations=10)
    # leid
}


### test grid search method with a recursive partition of cluster?
### choose AMI instead of ARI for unevenly distributed area label
### how to pursue 
    #- a well-partitioned area with 
    #- as much as possible consistence and low computation overhead overhead 
### to label? 
### `exp2.([0.1:0.2:10...])/1024`` for grid-search up to 1 





#’ calculate cross-specie entropy measure
#'
#' @param tb dataframe contain all information except 
#' @param kmclud k-means result
#' @param seed random seed
#' @param vln get tibble to plot vln plot
#' @return entropy calculation result
#' @noRd 
get_randomnized_entropy = function(tb, kmclud) {
    plan(multicore, workers=64)


    order = 1:nrow(kmclud)
    kmcludf = kmclud[sample(order),]

    cluster_entropy = function(cl) {
        rl = (cl+1)/sum(cl)
        -sum(rl*log2(rl))
    }
    leidf = future_sapply(1:ncol(kmcludf), function(g) {
        membership = kmcludf[,g]
        # df = as.data.frame(table(tb$specie, membership))
        df_mat = table(tb$specie, membership)
        df_mat = apply(df_mat, 2, cluster_entropy)
        df = tibble(entropy=df_mat, cluster=as.integer(names(df_mat)), reso=resos[g])
        # colnames(df) = c("cluster", "specie", "count")
        # df = df |> 
        #     mutate(cluster=as.integer(cluster)) |>
        #     group_by(cluster) |>
        #     summarise(entropy=cluster_entropy(count)) |>
        #     ungroup() |>
        #     mutate(reso=resos[g])
        mapping = setNames(df_mat, names(df_mat))
        mapping[kmcludf[,g]]
        # df = left_join(cbind(tb |> rownames_to_column("cell"), cluster=membership), df)
        # df

        # tb |> rownames_to_column("cell") |> select(cell, specie, chip, mainAreaID, ld, entropy, cluster, reso, specie_maid)
    })


    leidfal = tb |> rownames_to_column("cell") |> select(cell, specie, chip, mainAreaID, specie_maid)
    leidfal$cell_entropy = rowMeans(leidf)
    matl = lapply(c("mouse", "macaque", "marmoset"),
        function(sp) {
            leidfals = leidfal |> 
                    filter(specie==!!sp) |>
                    mutate(entropy=cell_entropy) |>
                    group_by(chip, mainAreaID) |>
                    summarise(entropy=exp(mean(cell_entropy))) |>
                    ungroup()
            mat = acast(leidfals, mainAreaID~chip, value.var="entropy")
            mat_na = setdiff(all_map, rownames(mat))
            mat = rbind(mat,
                matrix(data=NA, ncol=ncol(mat), nrow=length(mat_na), dimnames=list(mat_na, colnames(mat))))[uarea,]
            mat[uarea,]
    })
    return(list(vln=leidfal$cell_entropy, mat=Reduce(cbind, matl)))
}


# COLORS = ((244,  67,  54),
#           (233,  30,  99),
#           (156,  39, 176),
#           (103,  58, 183),
#           ( 63,  81, 181),
#           ( 33, 150, 243),
#           (  3, 169, 244),
#           (  0, 188, 212),
#           (  0, 150, 136),
#           ( 76, 175,  80),
#           (139, 195,  74),
#           (205, 220,  57),
#           (255, 235,  59),
#           (255, 193,   7),
#           (255, 152,   0),
#           (255,  87,  34),
#           (121,  85,  72),
#           (158, 158, 158),
#           ( 96, 125, 139))




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
