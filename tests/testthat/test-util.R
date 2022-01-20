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

cnmat <- subset(createCNmatrix(sim_data_bb$ascn), select = -c(chr, start, end, width))
cnmat_bp1 <- subset(createbreakpointmatrix(sim_data_bb$ascn)$bps, select = -c(loci))
cnmat_bp2 <- subset(createbreakpointmatrix(sim_data_bb$ascn, use_state = TRUE)$bps, select = -c(loci))
cnmat_bp3 <- subset(createbreakpointmatrix(sim_data_bb$ascn, use_state = TRUE, internalonly = TRUE)$bps, select = -c(loci))


test_that("Test matrix dimensions", {
  expect_equal(dim(cnmat)[2], sum(clones_dist))
  expect_equal(dim(cnmat_bp1)[2], sum(clones_dist))
  expect_equal(dim(cnmat_bp2)[2], sum(clones_dist))
  expect_equal(dim(cnmat_bp3)[2], sum(clones_dist))
  expect_equal(Mode(c(1, 2, 2, 2, 3, 4, 5)), 2)
})

wide_haps <- widen_haplotypebins(sim_data_bb$haplotypes, binsize = 5e6)
wide_cn <- widen_bins(sim_data_bb$CNbins, binsize = 5e6)
test_that("Test widen bins functions", {
  expect_equal(wide_haps$end[1] - wide_haps$start[1], 5e6 - 1)
  expect_equal(wide_cn$end[1] - wide_cn$start[1], 5e6 - 1)
})

segs <- create_segments(sim_data_bb$CNbins)
test_that("Test create segments", {
  expect_equal(dim(segs)[1], 24 * sum(clones_dist))
})

segs_cn <- create_segments(consensuscopynumber(CNbins))
segs_filt <- filter_segments(segs_cn, binwidth = 10e6)
test_that("Test filtering segments", {
  expect_lt(dim(segs_filt)[1], dim(segs_cn)[1])
  expect_true(segs_filt %>% dplyr::mutate(w = end - start) %>% dplyr::pull(w) %>% min() > 10e6)
})

bins <- segments_to_bins(segs_cn, binsize = 0.5e6)
bin_cell_id <- paste(bins$cell_id, bins$chr, bins$start, bins$end, sep = "_")
test_that("Test segments to bins", {
  expect_lt(dim(segs_cn)[1], dim(bins)[1])
  expect_equal(median((bins$end - bins$start) + 1), 0.5e6)
  expect_true(all(!duplicated(bin_cell_id))) #no duplicates
})
