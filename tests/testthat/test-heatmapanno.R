sim_data_bb <- simulate_data_cohort(
  clone_num = c(20, 25, 25, 10),
  clonal_events = list(
    list("1" = c(2, 0), "5" = c(3, 1)),
    list("2" = c(6, 3), "3" = c(1, 0)),
    list("17" = c(3, 1), "8" = c(6, 2)),
    list("1" = c(2, 2), "9" = c(4, 1))
  ), # opposite LOH on chr 1
  loherror = 0.01,
  coverage = 100,
  rho = 0.02,
  likelihood = "betabinomial",
  nchr = 0
)

dfannot <- data.frame(cell_id = unique(sim_data_bb$ascn$cell_id)) 
dfannot$Event <- NULL
dfannot$Event[1:20] <- "1"
dfannot$Event[21:45] <- "2"
dfannot$Event[46:70] <- "3"
dfannot$Event[71:80] <- "4"
dfannot$Clone <- "0"
dfannot$Other <- NULL
dfannot$Other[1:40] <- "Other1"
dfannot$Other[41:80] <- "Other2"

hm1 <- plotHeatmap(sim_data_bb$ascn, 
            annotations = dfannot, 
            tree = NULL, 
            reorderclusters = TRUE, 
            plottree = FALSE)

cl <- umap_clustering(sim_data_bb$ascn, 
                     field = "state",
                     minPts = 5,
                     umapmetric = "euclidean")
dfannot2 <- dplyr::left_join(dfannot, 
                            cl$clustering[c("cell_id", "clone_id")], 
                            by = "cell_id")
hm2 <- plotHeatmap(sim_data_bb$ascn, 
            annotations = dfannot2, 
            clusters = cl$clustering,
            tree = cl$tree, 
            reorderclusters = TRUE, 
            plottree = TRUE)

test_that("Test returns plot object", {
  expect_true(typeof(hm1) == "S4")
  expect_true(typeof(hm2) == "S4")
})
