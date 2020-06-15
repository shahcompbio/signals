library(ggplot2)
loherror <- 0.02
nclones <- 4
clones_dist <- c(20, 15, 35, 10)
sim_data_bb <- simulate_data_cohort(clone_num = clones_dist,
                                    clonal_events = list(list("1" = c(3,0), "5" = c(3,1), "6" = c(3,1)),
                                                         list("2" = c(6,3), "3" = c(1,0), "6" = c(1,0)),
                                                         list("17" = c(3,1), "8" = c(6,2), "19" = c(5,1)),
                                                         list("4" = c(5,2), "9" = c(4,1)), "2" = c(3,1)), #opposite LOH on chr 1
                                    loherror = loherror,
                                    coverage = 30,
                                    rho = 0.02,
                                    likelihood = "betabinomial",
                                    nchr = 0)

cl_state <- umap_clustering(sim_data_bb$CNbins, field = "copy", minPts = 5)
cl_bps <- umap_clustering_breakpoints(sim_data_bb$CNbins, minPts = 5,
                                      internalonly = FALSE, use_state = TRUE, state_remove = 2, n_neighbors = 10)
plot1 <- plot_umap(cl_state$clustering)

test_that("Test umap clustering with bin states", {
  expect_equal(length(unique(cl_state$clustering$clone_id)), 4)
  expect_true(all(table(cl_state$clustering$clone_id) %in% clones_dist))
  expect_true(is.ggplot(plot1))
})

test_that("Test umap clustering with breakpoints matrix", {
  expect_equal(length(unique(cl_bps$clustering$clone_id)), 4)
  expect_true(all(table(cl_bps$clustering$clone_id) %in% clones_dist))
})




