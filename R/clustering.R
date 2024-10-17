#' Perform UMAP dimensionality reduction and HDBSCAN clustering on copy number data
#'
#' This function takes copy number data, performs UMAP dimensionality reduction,
#' and then applies HDBSCAN clustering to identify cell populations. It can handle
#' both standard copy number data and haplotype-specific copy number (HSCN) data.
#'
#' @param CNbins A data frame containing copy number data. Must include columns
#'   for 'cell_id' and the specified `field`.
#' @param n_neighbors Integer. The number of neighbors to consider in UMAP. Default is 10.
#' @param min_dist Numeric. The minimum distance between points in UMAP. Default is 0.1.
#' @param minPts Integer. The minimum number of points to form a cluster in HDBSCAN. Default is 30.
#' @param seed Integer or NULL. Random seed for reproducibility. Default is NULL.
#' @param field Character. The column name in `CNbins` to use for copy number values. Default is "copy".
#' @param umapmetric Character. The distance metric to use in UMAP. Default is "correlation".
#' @param hscn Logical. Whether to use haplotype-specific copy number data. Default is FALSE.
#' @param pca Integer or NULL. Number of principal components to use in UMAP.  If NULL, pca not used, this is the default.
#'
#'
#' @return A list containing:
#'   \item{clustering}{A data frame with UMAP coordinates and cluster assignments for each cell.}
#'   \item{hdbscanresults}{The results of the HDBSCAN clustering.}
#'   \item{umapresults}{The results of the UMAP dimensionality reduction.}
#'   \item{tree}{A phylogenetic tree object representing the hierarchical structure of the clusters.}
#'
#' @details
#' The function performs the following steps:
#' 1. Creates a copy number matrix from the input data.
#' 2. Applies UMAP dimensionality reduction.
#' 3. Performs HDBSCAN clustering on the UMAP results.
#' 4. Generates a phylogenetic tree from the clustering results.
#'
#' If `hscn` is TRUE, the function expects columns 'copy' and 'BAF' in `CNbins`,
#' and creates separate matrices for A and B alleles.
#'
#' The function automatically adjusts `n_neighbors` if there are too few cells.
#' If UMAP fails, it attempts to rerun with small jitter added to the data points.
#' The function will reduce `minPts` if only one cluster is initially found.
#'
#' @export
umap_clustering <- function(CNbins,
                            n_neighbors = 10,
                            min_dist = 0.1,
                            minPts = 30,
                            seed = NULL,
                            field = "copy",
                            umapmetric = "correlation",
                            hscn = FALSE,
                            pca = NULL) {
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
  if (nrow(cnmatrix) > 500 & is.null(seed) & !is.null(pca)) {
    pca <- min(pca, ncol(cnmatrix))
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
                                fast_sgd = fast_sgd,
                                pca_method = "svdr")
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
                                fast_sgd = fast_sgd,
                                pca_method = "svdr")
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
