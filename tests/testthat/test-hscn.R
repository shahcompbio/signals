loherror <- 0.02
sim_data_bb <- simulate_data_cohort(clone_num = c(20, 25, 25, 10),
                                    clonal_events = list(list("1" = c(2,0), "5" = c(3,1)),
                                                         list("2" = c(6,3), "3" = c(1,0)),
                                                         list("17" = c(3,1), "8" = c(6,2)),
                                                         list("1" = c(2,2), "9" = c(4,1))), #opposite LOH on chr 1
                                    loherror = loherror,
                                    coverage = 30,
                                    rho = 0.02,
                                    likelihood = "betabinomial",
                                    nchr = 0)
results_bb <- callHaplotypeSpecificCN(sim_data_bb$CNbins, sim_data_bb$haplotypes, likelihood = "betabinomial")
f <- table(results_bb$data[results_bb$data$chr == "1",]$state_AS_phased)
f <- f / sum(f)

results_df <- orderdf(results_bb$data)
truth_df <- orderdf(sim_data_bb$ascn)
nbins <- length(truth_df$cell_id)

test_that("Test haplotype specific copy number inference (beta-binomial)", {
  expect_gt(sum(results_df$state_AS_phased == truth_df$state_AS_phased) / nbins, 0.99) #test accuraccy is > 99%
  expect_equal(round(results_bb$loherror, 2), loherror)
  expect_gt(results_bb$likelihood$taronesZ, 5)
  expect_true(f[["0|2"]] == 0.125)
  expect_true(f[["2|0"]] == 0.25)
})
