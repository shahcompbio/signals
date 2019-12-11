#' @export
umap_clustering <- function(CNbins,
                            n_neighbors = 30,
                            min_dist = 0.0,
                            minPts = 30,
                            seed = 1,
                            field = "state"){

  message("Creating CN matrix...")
  cnmatrix <- createCNmatrix(CNbins, na.rm = TRUE, field = field)
  cnmatrix <- subset(cnmatrix, select = -c(chr, start, end, width))
  cnmatrix <- t(cnmatrix)

  set.seed(seed)
  message('Calculating UMAP dimensionality reduction...')
  umapresults <- uwot::umap(cnmatrix,
                      n_neighbors = n_neighbors,
                      n_components = 2,
                      min_dist = min_dist,
                      ret_model = TRUE,
                      ret_nn = TRUE)

  dfumap <- data.frame(umap1 = umapresults$embedding[,1],
                       umap2 = umapresults$embedding[,2],
                       cell_id = row.names(cnmatrix))
  rownames(dfumap) <- row.names(cnmatrix)

  message('Clustering cells using hdbscan...')
  hdbscanresults <- dbscan::hdbscan(dfumap[,1:2], minPts = minPts,
                                    gen_simplified_tree = TRUE)
  clusterids <- hdbscanresults$cluster
  clusterids[clusterids == 0] <- 702
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  dfumap$clone_id <- LETTERS702[clusterids]
  dfumap <- dfumap %>%
    dplyr::mutate(clone_id = ifelse(clone_id == "ZZ", "0", clone_id))

  tree <- ape::as.phylo(hdbscanresults$hc, use.labels = TRUE)
  tree$tip.label <- row.names(cnmatrix)[as.numeric(tree$tip.label)]

  message(paste0("Identified ", length(unique(dfumap$clone_id)), " clusters"))
  message("Distribution of clusters:")
  print(table(dfumap$clone_id))

  return(list(clustering = dfumap,
              hdbscanresults = hdbscanresults,
              umapresults = umapresults,
              tree = tree))
}
