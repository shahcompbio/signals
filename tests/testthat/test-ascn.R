loherror <- 0.02
sim_data <- simulate_data_cohort(
  clone_num = c(5, 5, 5),
  clonal_events = list(
    list("1" = c(2, 0), "5" = c(3, 1)),
    list("2" = c(6, 3), "3" = c(1, 0)),
    list("1" = c(3, 1), "8" = c(6, 2))
  ),
  loherror = loherror,
  coverage = 30,
  nchr = 0
)
results <- callAlleleSpecificCN(sim_data$CNbins, sim_data$haplotypes)
results_2 <- callAlleleSpecificCN(sim_data$CNbins, sim_data$haplotypes, likelihood = "betabinomial")

print(results)

results_df <- orderdf(results$data)
truth_df <- orderdf(sim_data$ascn)
nbins <- length(truth_df$cell_id)

test_that("Test allele specific copy number inference", {
  expect_gt(sum(results_df$state_AS_phased == truth_df$state_AS_phased) / nbins, 0.99) # test accuraccy is > 99%
  expect_equal(round(results$loherror, 2), loherror)
  expect_lt(results_2$likelihood$taronesZ, 5)
})


loherror <- 0.02
sim_data_bb <- simulate_data_cohort(
  clone_num = c(5, 5, 5),
  clonal_events = list(
    list("1" = c(2, 0), "5" = c(3, 1)),
    list("2" = c(6, 3), "3" = c(1, 0)),
    list("1" = c(3, 1), "8" = c(6, 2))
  ),
  loherror = loherror,
  coverage = 30,
  rho = 0.02,
  likelihood = "betabinomial",
  nchr = 1
)
results_bb <- callAlleleSpecificCN(sim_data_bb$CNbins, sim_data_bb$haplotypes, likelihood = "betabinomial")

results_bb_df <- orderdf(results_bb$data)
truth_bb_df <- orderdf(sim_data_bb$ascn)
nbins <- length(truth_bb_df$cell_id)

test_that("Test allele specific copy number inference (beta-binomial)", {
  expect_gt(sum(results_bb_df$state_AS_phased == truth_bb_df$state_AS_phased) / nbins, 0.99) # test accuraccy is > 99%
  expect_equal(round(results_bb$loherror, 2), loherror)
  expect_gt(results_bb$likelihood$taronesZ, 5)
})


test_that("alleleHMM tolerates homozygous-deletion bins without corrupting the rest of the chromosome", {
  # A homdel bin (total CN 0) must not produce all-(-Inf) emissions, which would
  # propagate through the Viterbi recursion and force every downstream bin to 0.
  # Bins 4-5 are identical to bins 1-2, so they must decode the same.
  n         <- c(40, 40,  0, 40, 40)
  x         <- c(20, 20,  0, 20, 20)
  binstates <- c( 2,  2,  0,  2,  2)
  minor_cn  <- 0:4

  res <- alleleHMM(n, x, binstates, minor_cn, likelihood = "binomial")

  # no bin should have all non-finite emissions
  all_inf <- apply(res$l, 1, function(row) all(!is.finite(row)))
  expect_false(any(all_inf))
  # downstream bins keep their balanced call rather than collapsing to 0
  expect_equal(res$minorcn, c(1, 1, 0, 1, 1))
})
