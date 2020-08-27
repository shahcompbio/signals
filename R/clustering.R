#' @export
umap_clustering <- function(CNbins,
                            n_neighbors = 10,
                            min_dist = 0.1,
                            minPts = 30,
                            seed = 1,
                            field = "state",
                            umapmetric = "euclidean"){

  if(length(unique(CNbins$cell_id)) < n_neighbors) {
    n_neighbors <- length(unique(CNbins$cell_id)) - 1
  }

  minPts <- max(minPts, 2)

  message("Creating CN matrix...")
  cnmatrix <- createCNmatrix(CNbins, fillna = TRUE, field = field)
  cnmatrix <- subset(cnmatrix, select = -c(chr, start, end, width))
  cnmatrix <- t(cnmatrix)
  cnmatrix[!is.finite(cnmatrix)] <- 0 # remove non finite values

  set.seed(seed)
  message('Calculating UMAP dimensionality reduction...')
  if (nrow(cnmatrix) > 500){
    pca <- 50
    fast_sgd <- TRUE
  } else{
    pca <- NULL
    fast_sgd <- FALSE
  }
  umapresults <- uwot::umap(cnmatrix,
                      metric = umapmetric,
                      n_neighbors = n_neighbors,
                      n_components = 2,
                      min_dist = min_dist,
                      ret_model = TRUE,
                      ret_nn = TRUE,
                      pca = pca,
                      fast_sgd = fast_sgd)

  dfumap <- data.frame(umap1 = umapresults$embedding[,1],
                       umap2 = umapresults$embedding[,2],
                       cell_id = row.names(cnmatrix))
  rownames(dfumap) <- row.names(cnmatrix)

  message('Clustering cells using hdbscan...')
  gentree <- FALSE
  while(gentree == FALSE){
    hdbscanresults <- try(dbscan::hdbscan(dfumap[,1:2], minPts = minPts,
                                      gen_hdbscan_tree = FALSE,
                                      gen_simplified_tree = FALSE))
    if (class(hdbscanresults) == "try-error"){
      message("Only 1 cluster found, reducing minPts size by 10...")
      minPts <- round(minPts - 10)
      if (minPts < 0){
        message("Only 1 cluster can be found")
        gentree <- TRUE
      }
      message(paste0("Cluster size = ", minPts))
    } else{
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
  for (cl in sort(unique(dfumap$clone_id))){
    message(paste0("  Cluster ", cl, ": ", f[[cl]]))
  }

  if (length(unique(dfumap$clone_id)) == 1){
    dfumap$clone_id <- "A"
  }

  tree <- ape::as.phylo(hdbscanresults$hc, use.labels = TRUE)
  tree$tip.label <- row.names(cnmatrix)[as.numeric(tree$tip.label)]

  return(list(clustering = dfumap,
              hdbscanresults = hdbscanresults,
              umapresults = umapresults,
              tree = tree))
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
                            fixjitter = TRUE){

  if(length(unique(CNbins$cell_id)) < n_neighbors) {
    n_neighbors <- length(unique(CNbins$cell_id)) - 1
  }

  minPts <- max(minPts, 2)

  message("Creating breakpoint matrix...")
  segs <- schnapps::create_segments(CNbins, field = field)
  segs_matrix <- createbreakpointmatrix(segs,
                                        transpose = TRUE,
                                        internalonly = internalonly,
                                        use_state = use_state,
                                        state_remove = state_remove,
                                        fixjitter = fixjitter)
  segs_matrix <- segs_matrix$bps

  set.seed(seed)
  message('Calculating UMAP dimensionality reduction...')

  if (nrow(segs_matrix) > 500){
    fast_sgd <- TRUE
  } else{
    fast_sgd <- FALSE
  }

  umapresults <- uwot::umap(segs_matrix,
                            metric = "hamming",
                            n_neighbors = n_neighbors,
                            n_components = 2,
                            min_dist = min_dist,
                            ret_model = TRUE,
                            ret_nn = TRUE,
                            fast_sgd = fast_sgd)

  dfumap <- data.frame(umap1 = umapresults$embedding[,1],
                       umap2 = umapresults$embedding[,2],
                       cell_id = row.names(segs_matrix))
  rownames(dfumap) <- row.names(segs_matrix)
  message('Clustering cells using hdbscan...')
  gentree <- FALSE
  while(gentree == FALSE){
    hdbscanresults <- try(dbscan::hdbscan(dfumap[,1:2], minPts = minPts,
                                          gen_simplified_tree = FALSE))
    if (class(hdbscanresults) == "try-error"){
      message("Only 1 cluster found, reducing minPts size by 10...")
      minPts <- round(minPts - 10)
      message(paste0("Cluster size = ", minPts))
      if (minPts < 0){
        message("Only 1 cluster can be found")
        gentree <- TRUE
      }
    } else{
      gentree <- TRUE
    }
  }
  clusterids <- hdbscanresults$cluster
  clusterids[clusterids == 0] <- 702
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  dfumap$clone_id <- LETTERS702[clusterids]
  dfumap <- dfumap %>%
    dplyr::mutate(clone_id = ifelse(clone_id == "ZZ", "0", clone_id))

  if (length(unique(dfumap$clone_id)) == 1){
    dfumap$clone_id <- "A"
  }

  tree <- ape::as.phylo(hdbscanresults$hc, use.labels = TRUE)
  tree$tip.label <- row.names(segs_matrix)[as.numeric(tree$tip.label)]

  message(paste0("Identified ", length(unique(dfumap$clone_id)), " clusters"))
  message("Distribution of clusters:")
  f <- table(dfumap$clone_id)
  for (cl in sort(unique(dfumap$clone_id))){
    message(paste0("  Cluster ", cl, ": ", f[[cl]]))
  }

  return(list(clustering = dfumap,
              hdbscanresults = hdbscanresults,
              umapresults = umapresults,
              tree = tree))
}
