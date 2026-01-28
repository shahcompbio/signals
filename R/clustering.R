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
      warning("UMAP calculation failed: ", e$message, "\n",
              "Retrying with small jitter added to data. Results may differ slightly.",
              call. = FALSE)

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

#' Perform PCA + kNN graph + Leiden clustering on copy number data
#'
#' This function performs dimensionality reduction via PCA, constructs a k-nearest
#' neighbor graph, and applies the Leiden community detection algorithm to identify
#' cell populations. It can handle both standard copy number data and haplotype-specific
#' copy number (HSCN) data.
#'
#' Inspired by community detection approaches developed by Sohrab Salehi.
#' TODO: Add reference to Salehi et al. paper on community detection in single-cell genomics.
#'
#' @param CNbins A data frame containing copy number data. Must include columns
#'   for 'cell_id' and the specified `field`.
#' @param field Character. The column name in `CNbins` to use for copy number values. Default is "copy".
#' @param n_pcs Integer. The number of principal components to compute. Default is 50.
#' @param k Integer. The number of nearest neighbors for graph construction. Default is 15.
#' @param resolution Numeric. Resolution parameter for Leiden algorithm (higher = more clusters). Default is 0.7.
#' @param z_clip Numeric. Maximum absolute z-score for clipping scaled data. Default is 10.
#' @param seed Integer or NULL. Random seed for reproducibility. Default is NULL.
#' @param hscn Logical. Whether to use haplotype-specific copy number data. Default is FALSE.
#' @param objective_function Character. Leiden objective function: "modularity" or "CPM". Default is "modularity".
#' @param tree_type Character. Type of phylogenetic tree to generate: "centroid" (flat clusters) or "cell" (hierarchical within clusters). Default is "centroid".
#'
#' @return A list containing:
#'   \item{clustering}{A data frame with cell_id and clone_id (cluster assignments).}
#'   \item{leiden_results}{The igraph communities object from Leiden clustering.}
#'   \item{pca_results}{The prcomp object from PCA.}
#'   \item{tree}{A phylogenetic tree object (cluster-level or cell-level based on tree_type).}
#'
#' @details
#' The function performs the following steps:
#' 1. Creates a copy number matrix from the input data.
#' 2. Applies z-score standardization with clipping to handle outliers.
#' 3. Performs PCA dimensionality reduction.
#' 4. Constructs a symmetric k-nearest neighbor graph in PCA space.
#' 5. Applies Leiden community detection algorithm.
#' 6. Generates a phylogenetic tree (either cluster centroids or full cell hierarchy).
#'
#' If `hscn` is TRUE, the function expects columns 'copy' and 'BAF' in `CNbins`,
#' and creates separate matrices for A and B alleles.
#'
#' The function automatically adjusts `k` if there are too few cells.
#' Unlike HDBSCAN (used in umap_clustering), Leiden produces flat cluster assignments,
#' so tree generation uses hierarchical clustering on cluster centroids to create a backbone,
#' then grafts cell subtrees onto it. Both tree types preserve clone blocks in the tree structure:
#' - `tree_type = "centroid"`: Cells within each cluster form a flat star/polytomy (no within-cluster hierarchy)
#' - `tree_type = "cell"`: Cells within each cluster are hierarchically organized via hclust
#'
#' @export
leiden_clustering <- function(CNbins,
                               field = "copy",
                               n_pcs = 50,
                               k = 15,
                               resolution = 0.7,
                               z_clip = 10,
                               seed = NULL,
                               hscn = FALSE,
                               objective_function = "modularity",
                               tree_type = "centroid") {
  # Check for igraph package
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package \"igraph\" needed for this function. Please install it.", call. = FALSE)
  }

  # Validate objective function
  if (!objective_function %in% c("modularity", "CPM")) {
    stop("objective_function must be 'modularity' or 'CPM'", call. = FALSE)
  }

  # Validate tree type
  if (!tree_type %in% c("centroid", "cell")) {
    stop("tree_type must be 'centroid' or 'cell'", call. = FALSE)
  }

  # Check field exists
  if (!field %in% names(CNbins)) {
    stop(paste0("Field '", field, "' not found in CNbins"), call. = FALSE)
  }

  n_cells <- length(unique(CNbins$cell_id))

  # Adjust k if necessary
  if (n_cells <= k) {
    k <- max(1, n_cells - 1)
    message(paste0("Adjusted k to ", k, " (n_cells - 1)"))
  }

  message("Creating CN matrix...")
  if (hscn) {
    CNbins$A <- CNbins$copy * CNbins$BAF
    CNbins$B <- CNbins$copy * (1 - CNbins$BAF)

    cnmatrixA <- createCNmatrix(CNbins, fillna = TRUE, field = "A")
    cnmatrixA <- subset(cnmatrixA, select = -c(chr, start, end, width))
    cnmatrixA <- t(cnmatrixA)
    cnmatrixA[!is.finite(cnmatrixA)] <- 0

    cnmatrixB <- createCNmatrix(CNbins, fillna = TRUE, field = "B")
    cnmatrixB <- subset(cnmatrixB, select = -c(chr, start, end, width))
    cnmatrixB <- t(cnmatrixB)
    cnmatrixB[!is.finite(cnmatrixB)] <- 0

    colnames(cnmatrixA) <- paste0(colnames(cnmatrixA), "_A")
    colnames(cnmatrixB) <- paste0(colnames(cnmatrixB), "_B")

    cnmatrix <- cbind(cnmatrixA, cnmatrixB)
  } else {
    cnmatrix <- createCNmatrix(CNbins, fillna = TRUE, field = field)
    cnmatrix <- subset(cnmatrix, select = -c(chr, start, end, width))
    cnmatrix <- t(cnmatrix)
    cnmatrix[!is.finite(cnmatrix)] <- 0
  }

  message("Scaling and clipping data...")
  cnmatrix_scaled <- scale(cnmatrix, center = TRUE, scale = TRUE)
  cnmatrix_scaled[!is.finite(cnmatrix_scaled)] <- 0
  cnmatrix_scaled <- pmax(pmin(cnmatrix_scaled, z_clip), -z_clip)

  message("Performing PCA dimensionality reduction...")
  n_pcs_actual <- min(n_pcs, nrow(cnmatrix_scaled) - 1, ncol(cnmatrix_scaled))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  pca_results <- stats::prcomp(cnmatrix_scaled, center = FALSE, scale. = FALSE, rank. = n_pcs_actual)
  pca_coords <- pca_results$x

  message("Building k-nearest neighbor graph...")
  # Compute Euclidean distances in PCA space
  dist_matrix <- as.matrix(dist(pca_coords, method = "euclidean"))

  # Build symmetric adjacency matrix from k nearest neighbors
  adj_matrix <- matrix(0, nrow = n_cells, ncol = n_cells)
  rownames(adj_matrix) <- rownames(pca_coords)
  colnames(adj_matrix) <- rownames(pca_coords)

  for (i in 1:n_cells) {
    # Find k nearest neighbors (excluding self)
    neighbors <- order(dist_matrix[i, ])[2:(k + 1)]
    adj_matrix[i, neighbors] <- 1
  }

  # Make symmetric
  adj_matrix <- adj_matrix + t(adj_matrix)
  adj_matrix[adj_matrix > 0] <- 1

  # Create igraph object
  knn_graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)

  message("Performing Leiden clustering...")
  if (!is.null(seed)) {
    set.seed(seed)
  }
  leiden_results <- igraph::cluster_leiden(
    knn_graph,
    objective_function = objective_function,
    resolution = resolution,
    n_iterations = 10
  )

  # Get cluster memberships
  membership <- igraph::membership(leiden_results)

  # Convert to letter labels
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  clone_ids <- LETTERS702[membership]

  # Create clustering dataframe
  clustering_df <- data.frame(
    cell_id = rownames(pca_coords),
    clone_id = clone_ids,
    stringsAsFactors = FALSE
  )
  rownames(clustering_df) <- clustering_df$cell_id

  message(paste0("Identified ", length(unique(clustering_df$clone_id)), " clusters"))
  message("Distribution of clusters:")
  f <- table(clustering_df$clone_id)
  for (cl in sort(unique(clustering_df$clone_id))) {
    message(paste0("  Cluster ", cl, ": ", f[[cl]]))
  }

  # Generate phylogenetic tree
  message(paste0("Generating ", tree_type, " phylogenetic tree..."))

  # Check for single cluster case
  cluster_labels <- sort(unique(clustering_df$clone_id))

  if (length(cluster_labels) == 1) {
    # Only one cluster - no need for backbone tree, just build tree of all cells
    message("Only one cluster found, creating simple tree of all cells...")
    all_cells <- clustering_df$cell_id

    if (tree_type == "centroid") {
      # Flat tree - all cells at same level
      if (length(all_cells) == 1) {
        tree <- ape::read.tree(text = paste0("(", all_cells, ");"))
      } else {
        newick_str <- paste0("(", paste(all_cells, collapse = ","), ");")
        tree <- ape::read.tree(text = newick_str)
      }
    } else {
      # tree_type == "cell" - hierarchical clustering of all cells
      if (length(all_cells) == 1) {
        tree <- ape::read.tree(text = paste0("(", all_cells, ");"))
      } else {
        hc <- stats::hclust(stats::dist(pca_coords), method = "average")
        tree <- ape::as.phylo(hc)
        tree$tip.label <- rownames(pca_coords)
      }
    }
  } else if (tree_type == "centroid") {
    # Multiple clusters - centroid tree with flat subtrees
    # Compute cluster centroids in PCA space
    centroids <- matrix(0, nrow = length(cluster_labels), ncol = ncol(pca_coords))
    rownames(centroids) <- cluster_labels

    for (i in seq_along(cluster_labels)) {
      cl <- cluster_labels[i]
      cells_in_cluster <- clustering_df$cell_id[clustering_df$clone_id == cl]
      centroids[i, ] <- colMeans(pca_coords[cells_in_cluster, , drop = FALSE])
    }

    # Build backbone tree from centroids
    hc_backbone <- stats::hclust(stats::dist(centroids), method = "average")
    backbone_tree <- ape::as.phylo(hc_backbone)
    backbone_tree$tip.label <- cluster_labels

    # Create flat/star subtrees for each cluster (no hierarchy within cluster)
    cell_subtrees <- list()
    for (cl in cluster_labels) {
      cells_in_cluster <- clustering_df$cell_id[clustering_df$clone_id == cl]
      if (length(cells_in_cluster) == 1) {
        # Single cell - create trivial tree
        subtree <- ape::read.tree(text = paste0("(", cells_in_cluster, ");"))
      } else {
        # Multiple cells - create star/polytomy (all cells at same level)
        newick_str <- paste0("(", paste(cells_in_cluster, collapse = ","), ");")
        subtree <- ape::read.tree(text = newick_str)
      }
      cell_subtrees[[cl]] <- subtree
    }

    # Graft flat subtrees onto backbone
    tree <- backbone_tree
    tip_order <- order(match(tree$tip.label, cluster_labels), decreasing = TRUE)

    for (idx in tip_order) {
      cl <- tree$tip.label[idx]
      subtree <- cell_subtrees[[cl]]
      tree <- ape::bind.tree(tree, subtree, where = idx)
    }
  } else {
    # tree_type == "cell" with multiple clusters
    # Graft cell subtrees onto centroid backbone to preserve cluster blocks

    # Step 1: Build centroid backbone tree
    centroids <- matrix(0, nrow = length(cluster_labels), ncol = ncol(pca_coords))
    rownames(centroids) <- cluster_labels

    for (i in seq_along(cluster_labels)) {
      cl <- cluster_labels[i]
      cells_in_cluster <- clustering_df$cell_id[clustering_df$clone_id == cl]
      centroids[i, ] <- colMeans(pca_coords[cells_in_cluster, , drop = FALSE])
    }

    hc_backbone <- stats::hclust(stats::dist(centroids), method = "average")
    backbone_tree <- ape::as.phylo(hc_backbone)
    backbone_tree$tip.label <- cluster_labels

    # Step 2: Build cell subtrees for each cluster
    cell_subtrees <- list()
    for (cl in cluster_labels) {
      cells_in_cluster <- clustering_df$cell_id[clustering_df$clone_id == cl]
      if (length(cells_in_cluster) == 1) {
        # Single cell - create trivial tree with one tip
        subtree <- ape::read.tree(text = paste0("(", cells_in_cluster, ");"))
      } else {
        # Multiple cells - hierarchical clustering within cluster
        cluster_pca <- pca_coords[cells_in_cluster, , drop = FALSE]
        hc_sub <- stats::hclust(stats::dist(cluster_pca), method = "average")
        subtree <- ape::as.phylo(hc_sub)
        subtree$tip.label <- cells_in_cluster
      }
      cell_subtrees[[cl]] <- subtree
    }

    # Step 3: Graft subtrees onto backbone
    # Process in reverse order of tip indices to avoid index shifting issues
    tree <- backbone_tree
    tip_order <- order(match(tree$tip.label, cluster_labels), decreasing = TRUE)

    for (idx in tip_order) {
      cl <- tree$tip.label[idx]
      subtree <- cell_subtrees[[cl]]
      # bind.tree replaces the tip at position 'where' with the subtree
      tree <- ape::bind.tree(tree, subtree, where = idx)
    }
  }

  return(list(
    clustering = clustering_df,
    leiden_results = leiden_results,
    pca_results = pca_results,
    tree = tree
  ))
}

#' Hierarchical clustering (deprecated)
#'
#' This function is deprecated and not implemented. Use [umap_clustering()] instead.
#'
#' @return NULL (with deprecation warning)
#' @export
hc_clustering <- function() {
  .Deprecated("umap_clustering", msg = "hc_clustering() is not implemented. Use umap_clustering() instead.")
  invisible(NULL)
}
