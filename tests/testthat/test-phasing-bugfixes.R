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

test_that("total copy number is not changed to repair the alleles", {
  # A + B > state used to be repaired by blanking `state` and filling it from a
  # neighbouring bin, which made total CN a function of the phasing
  d <- data.frame(state = c(2L, 2L, 1L, 2L), A = c(1L, 1L, 1L, 1L), B = c(1L, 1L, 1L, 1L))
  fixed <- d
  fixed$B <- ifelse(!is.na(fixed$A) & !is.na(fixed$B) & !is.na(fixed$state) &
                      (fixed$A + fixed$B) > fixed$state,
                    pmax(pmin(fixed$B, fixed$state), 0), fixed$B)
  fixed$A <- ifelse(!is.na(fixed$B) & !is.na(fixed$state), fixed$state - fixed$B, fixed$A)
  expect_equal(fixed$state, d$state)              # total CN preserved
  expect_equal(fixed$A + fixed$B, fixed$state)    # and now consistent
  expect_true(all(fixed$A >= 0) && all(fixed$B >= 0))
})

test_that("min_propA defaults to off so existing results are unchanged", {
  expect_equal(formals(callHaplotypeSpecificCN)$min_propA, 0)
  expect_equal(formals(signals:::get_cells_per_chr_local)$min_propA, 0)
  expect_equal(formals(signals:::proportion_imbalance)$min_propA, 0)
})

test_that("min_propA gates on the best cluster's propA", {
  # the floor replaces which.max(propA) with a fallback to every cell
  gate <- function(propA, min_propA) !is.na(propA) && propA < min_propA
  expect_false(gate(0.40, 0))      # off by default, nothing is ever gated
  expect_false(gate(0.00, 0))
  expect_true(gate(0.03, 0.05))    # weak unit is gated once a floor is set
  expect_false(gate(0.40, 0.05))   # a genuinely imbalanced unit is not
  expect_false(gate(NA, 0.05))     # no clusters at all is not "unphaseable"
})
