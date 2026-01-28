library(ggplot2)
loherror <- 0.02
nclones <- 4
clones_dist <- c(20, 15, 35, 10)
sim_data_bb <- simulate_data_cohort(
  clone_num = clones_dist,
  clonal_events = list(list("1" = c(3, 0), "5" = c(3, 1), "6" = c(3, 1)),
    list("2" = c(6, 3), "3" = c(1, 0), "6" = c(1, 0)),
    list("17" = c(3, 1), "8" = c(6, 2), "19" = c(5, 1)),
    list("4" = c(5, 2), "9" = c(4, 1)),
    "2" = c(3, 1)
  ), # opposite LOH on chr 1
  loherror = loherror,
  coverage = 30,
  rho = 0.02,
  likelihood = "betabinomial",
  nchr = 0
)

cl_state <- umap_clustering(sim_data_bb$CNbins, field = "copy", minPts = 5)
cl_bps <- umap_clustering_breakpoints(sim_data_bb$CNbins,
  minPts = 5,
  internalonly = FALSE, use_state = TRUE,
  state_remove = 2, n_neighbors = 10, fixjitter = FALSE
)
plot1 <- plot_umap(cl_state$clustering)

test_that("Test umap clustering with bin states", {
  expect_equal(length(unique(cl_state$clustering$clone_id)), 4)
  expect_true(all(table(cl_state$clustering$clone_id) %in% clones_dist))
  expect_true(ggplot2::is_ggplot(plot1))
})

test_that("Test umap clustering with breakpoints matrix", {
  expect_equal(length(unique(cl_bps$clustering$clone_id)), 4)
  expect_true(all(table(cl_bps$clustering$clone_id) %in% clones_dist))
})

test_that("leiden_clustering returns expected structure", {
  skip_if_not_installed("igraph")
  cl <- leiden_clustering(sim_data_bb$CNbins,
    field = "copy",
    resolution = 0.5,
    k = 10,
    seed = 42
  )
  expect_true(all(c("cell_id", "clone_id") %in% names(cl$clustering)))
  expect_s3_class(cl$pca_results, "prcomp")
  expect_true(inherits(cl$leiden_results, "communities"))
  expect_s3_class(cl$tree, "phylo")
  expect_true(length(unique(cl$clustering$clone_id)) > 1)
})

test_that("leiden_clustering handles tree_type parameter", {
  skip_if_not_installed("igraph")

  # Test centroid tree (flat clusters)
  cl_centroid <- leiden_clustering(sim_data_bb$CNbins,
    field = "copy",
    tree_type = "centroid",
    k = 10,
    seed = 42
  )
  expect_s3_class(cl_centroid$tree, "phylo")
  # Centroid tree should have cell IDs as tip labels (all cells, but flat within clusters)
  expect_true(all(cl_centroid$tree$tip.label %in% cl_centroid$clustering$cell_id))
  # All cells should be represented
  expect_equal(length(cl_centroid$tree$tip.label), nrow(cl_centroid$clustering))

  # Test cell tree (hierarchical clusters)
  cl_cell <- leiden_clustering(sim_data_bb$CNbins,
    field = "copy",
    tree_type = "cell",
    k = 10,
    seed = 42
  )
  expect_s3_class(cl_cell$tree, "phylo")
  # Cell tree should have cell IDs as tip labels
  expect_true(all(cl_cell$tree$tip.label %in% cl_cell$clustering$cell_id))
  # All cells should be represented
  expect_equal(length(cl_cell$tree$tip.label), nrow(cl_cell$clustering))
})

test_that("leiden_clustering works with hscn mode", {
  skip_if_not_installed("igraph")

  # Add BAF column for hscn mode
  sim_data_with_baf <- sim_data_bb$CNbins
  sim_data_with_baf$BAF <- 0.5

  cl_hscn <- leiden_clustering(sim_data_with_baf,
    field = "copy",
    hscn = TRUE,
    k = 10,
    seed = 42
  )
  expect_true(all(c("cell_id", "clone_id") %in% names(cl_hscn$clustering)))
  expect_s3_class(cl_hscn$pca_results, "prcomp")
})

test_that("cell tree preserves clone blocks (cells from same cluster are contiguous)", {
  skip_if_not_installed("igraph")

  cl <- leiden_clustering(sim_data_bb$CNbins,
    field = "copy",
    tree_type = "cell",
    k = 10,
    seed = 42
  )

  # Get tree tip order (how cells appear in the tree from left to right)
  tree_order <- cl$tree$tip.label[cl$tree$edge[cl$tree$edge[, 2] <= length(cl$tree$tip.label), 2]]
  # ape stores tips in edge order, but we need plot order
  tree_order <- cl$tree$tip.label

  # Map each cell to its clone
  cell_to_clone <- setNames(cl$clustering$clone_id, cl$clustering$cell_id)
  clone_order <- cell_to_clone[tree_order]

  # Check that each clone forms a contiguous block

  # Use run-length encoding: if blocks are contiguous, each clone appears in exactly one run
  rle_clones <- rle(clone_order)
  clone_runs <- table(rle_clones$values)

  # Each clone should appear exactly once in the RLE (meaning all its cells are together)
  expect_true(all(clone_runs == 1),
    info = paste("Some clones are split in the tree. Clone run counts:",
                 paste(names(clone_runs), clone_runs, sep = "=", collapse = ", ")))
})

test_that("centroid tree preserves clone blocks (cells from same cluster are contiguous)", {
  skip_if_not_installed("igraph")

  cl <- leiden_clustering(sim_data_bb$CNbins,
    field = "copy",
    tree_type = "centroid",
    k = 10,
    seed = 42
  )

  # Get tree tip order
  tree_order <- cl$tree$tip.label

  # Map each cell to its clone
  cell_to_clone <- setNames(cl$clustering$clone_id, cl$clustering$cell_id)
  clone_order <- cell_to_clone[tree_order]

  # Check that each clone forms a contiguous block
  rle_clones <- rle(clone_order)
  clone_runs <- table(rle_clones$values)

  # Each clone should appear exactly once in the RLE (meaning all its cells are together)
  expect_true(all(clone_runs == 1),
    info = paste("Some clones are split in the centroid tree. Clone run counts:",
                 paste(names(clone_runs), clone_runs, sep = "=", collapse = ", ")))
})
