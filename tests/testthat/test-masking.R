
loherror <- 0.02
sim_data_bb <- simulate_data_cohort(
  clone_num = c(20, 25, 25, 10),
  clonal_events = list(
    list("1" = c(2, 0), "5" = c(3, 1)),
    list("2" = c(6, 3), "3" = c(1, 0)),
    list("17" = c(3, 1), "8" = c(6, 2)),
    list("1" = c(2, 2), "9" = c(4, 1))
  ), # opposite LOH on chr 1
  loherror = loherror,
  coverage = 100,
  rho = 0.02,
  likelihood = "betabinomial",
  nchr = 0
)
results_bb_nonmask <- callHaplotypeSpecificCN(sim_data_bb$CNbins, 
                                      sim_data_bb$haplotypes, 
                                      likelihood = "betabinomial")

maskedbins <- sim_data_bb$CNbins %>% 
  dplyr::select(chr, start, end) %>% 
  head(10)

results_bb_mask <- callHaplotypeSpecificCN(sim_data_bb$CNbins, 
                                              sim_data_bb$haplotypes, 
                                              likelihood = "betabinomial",
                                              maskedbins = maskedbins)
cell1 <- results_bb_mask$data$cell_id[1]

cell_mask <- dplyr::filter(results_bb_mask$data, cell_id == cell1)
cell_unmask <- dplyr::filter(results_bb_nonmask$data, cell_id == cell1)

test_that("Test masking 10 bins from chromosome 1", {
  expect_equal(nrow(cell_mask), nrow(cell_unmask))
  expect_equal(nrow(na.omit(cell_mask)), nrow(cell_unmask) - 10)
  expect_equal(sum(is.na(cell_mask$BAF)), 10)
  expect_equal(sum(is.na(cell_unmask$BAF)), 0)
})
