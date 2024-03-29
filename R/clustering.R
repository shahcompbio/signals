#' @export
umap_clustering <- function(CNbins,
                            n_neighbors = 10,
                            min_dist = 0.1,
                            minPts = 30,
                            seed = NULL,
                            field = "copy",
                            umapmetric = "correlation",
                            hscn = FALSE) {
  if (length(unique(CNbins$cell_id)) < n_neighbors) {
    n_neighbors <- length(unique(CNbins$cell_id)) - 1
  }

  minPts <- max(minPts, 2)

  message("Creating CN matrix...")
  if (hscn){
    CNbins$A <- CNbins$copy * CNbins$BAF
    CNbins$B <- CNbins$copy * (1 - CNbins$BAF)
    
    cnmatrixA <- createCNmatrix(CNbins, fillna = TRUE, field = "A")
    cnmatrixA <- subset(cnmatrixA, select = -c(chr, start, end, width))
    cnmatrixA <- t(cnmatrixA)
    cnmatrixA[!is.finite(cnmatrixA)] <- 0 # remove non finite values
    
    cnmatrixB <- createCNmatrix(CNbins, fillna = TRUE, field = "B")
    cnmatrixB <- subset(cnmatrixB, select = -c(chr, start, end, width))
    cnmatrixB <- t(cnmatrixB)
    cnmatrixB[!is.finite(cnmatrixB)] <- 0 # remove non finite values
    
    colnames(cnmatrixA) <- paste0(colnames(cnmatrixA), "_A")
    colnames(cnmatrixB) <- paste0(colnames(cnmatrixB), "_B")
    
    cnmatrix <- cbind(cnmatrixA, cnmatrixB)
    
  } else{
    cnmatrix <- createCNmatrix(CNbins, fillna = TRUE, field = field)
    cnmatrix <- subset(cnmatrix, select = -c(chr, start, end, width))
    cnmatrix <- t(cnmatrix)
    cnmatrix[!is.finite(cnmatrix)] <- 0 # remove non finite values
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }
  message("Calculating UMAP dimensionality reduction...")
  if (nrow(cnmatrix) > 500 & is.null(seed)) {
    pca <- min(50, ncol(cnmatrix))
    fast_sgd <- TRUE
  } else {
    pca <- NULL
    fast_sgd <- FALSE
  }
  # umapresults <- uwot::umap(cnmatrix,
  #   metric = umapmetric,
  #   n_neighbors = n_neighbors,
  #   n_components = 2,
  #   min_dist = min_dist,
  #   ret_model = TRUE,
  #   ret_nn = TRUE,
  #   pca = pca,
  #   fast_sgd = fast_sgd
  # )

  #TODO find out why umap gives an error for some cases, seems to be a new bug
  umapresults <- tryCatch(
    {
      umapresults <- uwot::umap(cnmatrix,
                                metric = umapmetric,
                                n_neighbors = n_neighbors,
                                n_components = 2,
                                min_dist = min_dist,
                                ret_model = TRUE,
                                ret_nn = TRUE,
                                pca = pca,
                                fast_sgd = fast_sgd)
    },
    error = function(e) {
      # Handle error by rerunning UMAP with different parameters
      message("An error occurred in umap calculation: ", e$message)
      message("Rerunning UMAP after adding small jitter to data points...")
      
      mat <- cnmatrix + matrix(runif(nrow(cnmatrix) * ncol(cnmatrix),
                                     min=-0.005, max=0.005), 
                               nrow=nrow(cnmatrix), ncol=ncol(cnmatrix))
      
      umapresults <- uwot::umap(mat,
                                metric = umapmetric,
                                n_neighbors = n_neighbors,
                                n_components = 2,
                                min_dist = min_dist,
                                ret_model = TRUE,
                                ret_nn = TRUE,
                                pca = pca,
                                fast_sgd = fast_sgd)
    }
  )
  
  
  
  dfumap <- data.frame(
    umap1 = umapresults$embedding[, 1],
    umap2 = umapresults$embedding[, 2],
    cell_id = row.names(cnmatrix)
  )
  dfumap$umap1 <- unlist(lapply(dfumap$umap1, function(x) ifelse(!is.finite(x), 0.0, x))) # remove non finite values
  dfumap$umap2 <- unlist(lapply(dfumap$umap2, function(x) ifelse(!is.finite(x), 0.0, x))) # remove non finite values
  rownames(dfumap) <- row.names(cnmatrix)

  message("Clustering cells using hdbscan...")
  gentree <- FALSE
  while (gentree == FALSE) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    hdbscanresults <- try(dbscan::hdbscan(dfumap[, 1:2],
      minPts = minPts,
      gen_hdbscan_tree = FALSE,
      gen_simplified_tree = FALSE
    ))
    if (class(hdbscanresults) == "try-error") {
      message("Only 1 cluster found, reducing minPts size by 10...")
      minPts <- round(minPts - 10)
      if (minPts <= 0) {
        message("Only 1 cluster can be found")
        gentree <- TRUE
      }
      message(paste0("Cluster size = ", minPts))
    } else {
      gentree <- TRUE
    }
  }
  clusterids <- hdbscanresults$cluster
  clusterids[clusterids == 0] <- 702
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  dfumap$clone_id <- LETTERS702[clusterids]
  dfumap <- dfumap %>%
    dplyr::mutate(clone_id = ifelse(clone_id == "ZZ", "0", clone_id))

  message(paste0("Identified ", length(unique(dfumap$clone_id)), " clusters"))
  message("Distribution of clusters:")
  f <- table(dfumap$clone_id)
  for (cl in sort(unique(dfumap$clone_id))) {
    message(paste0("  Cluster ", cl, ": ", f[[cl]]))
  }

  if (length(unique(dfumap$clone_id)) == 1) {
    dfumap$clone_id <- "A"
  }

  tree <- ape::as.phylo(hdbscanresults$hc, use.labels = TRUE)
  tree$tip.label <- row.names(cnmatrix)[as.numeric(tree$tip.label)]

  return(list(
    clustering = dfumap,
    hdbscanresults = hdbscanresults,
    umapresults = umapresults,
    tree = tree
  ))
}

#' @export
umap_clustering_breakpoints <- function(CNbins,
                                        n_neighbors = 10,
                                        min_dist = 0.1,
                                        minPts = 30,
                                        seed = 1,
                                        field = "state",
                                        internalonly = TRUE,
                                        use_state = FALSE,
                                        state_remove = 2,
                                        fixjitter = TRUE) {
  if (length(unique(CNbins$cell_id)) < n_neighbors) {
    n_neighbors <- length(unique(CNbins$cell_id)) - 1
  }

  minPts <- max(minPts, 2)

  message("Creating breakpoint matrix...")
  segs <- create_segments(CNbins, field = field)
  segs_matrix <- createbreakpointmatrix(segs,
    transpose = TRUE,
    internalonly = internalonly,
    use_state = use_state,
    state_remove = state_remove,
    fixjitter = fixjitter
  )
  segs_matrix <- segs_matrix$bps

  if (!is.null(seed)) {
    set.seed(seed)
  }
  message("Calculating UMAP dimensionality reduction...")

  if (nrow(segs_matrix) > 500 & is.null(seed)) {
    fast_sgd <- TRUE
  } else {
    fast_sgd <- FALSE
  }

  umapresults <- uwot::umap(segs_matrix,
    metric = "hamming",
    n_neighbors = n_neighbors,
    n_components = 2,
    min_dist = min_dist,
    ret_model = TRUE,
    ret_nn = TRUE,
    fast_sgd = fast_sgd
  )

  dfumap <- data.frame(
    umap1 = umapresults$embedding[, 1],
    umap2 = umapresults$embedding[, 2],
    cell_id = row.names(segs_matrix)
  )
  rownames(dfumap) <- row.names(segs_matrix)
  message("Clustering cells using hdbscan...")
  gentree <- FALSE
  while (gentree == FALSE) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    hdbscanresults <- try(dbscan::hdbscan(dfumap[, 1:2],
      minPts = minPts,
      gen_simplified_tree = FALSE
    ))
    if (class(hdbscanresults) == "try-error") {
      message("Only 1 cluster found, reducing minPts size by 10...")
      minPts <- round(minPts - 10)
      message(paste0("Cluster size = ", minPts))
      if (minPts < 0) {
        message("Only 1 cluster can be found")
        gentree <- TRUE
      }
    } else {
      gentree <- TRUE
    }
  }
  clusterids <- hdbscanresults$cluster
  clusterids[clusterids == 0] <- 702
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  dfumap$clone_id <- LETTERS702[clusterids]
  dfumap <- dfumap %>%
    dplyr::mutate(clone_id = ifelse(clone_id == "ZZ", "0", clone_id))

  if (length(unique(dfumap$clone_id)) == 1) {
    dfumap$clone_id <- "A"
  }

  tree <- ape::as.phylo(hdbscanresults$hc, use.labels = TRUE)
  tree$tip.label <- row.names(segs_matrix)[as.numeric(tree$tip.label)]

  message(paste0("Identified ", length(unique(dfumap$clone_id)), " clusters"))
  message("Distribution of clusters:")
  f <- table(dfumap$clone_id)
  for (cl in sort(unique(dfumap$clone_id))) {
    message(paste0("  Cluster ", cl, ": ", f[[cl]]))
  }

  return(list(
    clustering = dfumap,
    hdbscanresults = hdbscanresults,
    umapresults = umapresults,
    tree = tree
  ))
}

#' @export
hc_clustering <- function() {

}
