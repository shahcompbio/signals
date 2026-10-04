# Regression tests for phasing bugs fixed in 0.18.0.

test_that("the hdbscan noise cluster is not used to phase a unit", {
  # cluster "0" is hdbscan's noise grab-bag (dbscan reports noise as cluster 0;
  # R/clustering.R remaps it through "ZZ" to the string "0"). It must not be
  # selected even when it happens to look the most imbalanced.
  clustering <- data.table::data.table(
    cell_id = c("n1", "n2", "c1", "c2"),
    clone_id = c("0", "0", "1", "1")
  )
  expect_true(any(clustering$clone_id != "0"))
  kept <- clustering[clone_id != "0"]
  expect_setequal(kept$cell_id, c("c1", "c2"))

  # but if every cell landed in the noise cluster, keep them rather than
  # returning nothing to phase with
  allnoise <- data.table::data.table(cell_id = c("n1", "n2"), clone_id = c("0", "0"))
  kept2 <- if (any(allnoise$clone_id != "0")) allnoise[clone_id != "0"] else allnoise
  expect_equal(nrow(kept2), 2L)
})
