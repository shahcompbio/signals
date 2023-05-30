
loherror <- 0.02
sim_data_bb1 <- simulate_data_cohort(
  clone_num = c(20),
  clonal_events = list(
    list("1" = c(2, 0))),
  loherror = loherror,
  coverage = 20,
  rho = 0.02,
  likelihood = "betabinomial",
  nchr = 0
)

sim_data_bb2 <- simulate_data_cohort(
  clone_num = c(20),
  clonal_events = list(
    list("1" = c(2, 1))),
  loherror = loherror,
  coverage = 20,
  rho = 0.02,
  likelihood = "betabinomial",
  nchr = 0
)

sim_data_bb3 <- simulate_data_cohort(
  clone_num = c(20),
  clonal_events = list(
    list("1" = c(2, 0))),
  loherror = loherror,
  coverage = 20,
  rho = 0.02,
  likelihood = "betabinomial",
  nchr = 0
)

sim_data_bb4 <- simulate_data_cohort(
  clone_num = c(20),
  clonal_events = list(
    list("1" = c(2, 1))),
  loherror = loherror,
  coverage = 20,
  rho = 0.02,
  likelihood = "betabinomial",
  nchr = 0
)

start_switch <- 80e6
end_switch <- 150e6

sim_data_bb2$CNbins$cell_id <- sim_data_bb1$CNbins$cell_id
sim_data_bb2$haplotypes$cell_id <- sim_data_bb1$haplotypes$cell_id
sim_data_bb3$CNbins$cell_id <- sim_data_bb4$CNbins$cell_id
sim_data_bb3$haplotypes$cell_id <- sim_data_bb4$haplotypes$cell_id

cndata <- sim_data_bb1$CNbins %>% 
  dplyr::filter(start < start_switch | start > end_switch)
cndata <- cndata %>% 
  dplyr::bind_rows(sim_data_bb2$CNbins %>% 
  dplyr::filter(start > start_switch & start <= end_switch))
cndata <- cndata %>% 
  dplyr::bind_rows(sim_data_bb4$CNbins %>% 
                     dplyr::filter(start < start_switch | start > end_switch))
cndata <- cndata %>% 
  dplyr::bind_rows(sim_data_bb3$CNbins %>% 
                     dplyr::filter(start > start_switch & start <= end_switch))

hapsdata <- sim_data_bb1$haplotypes %>% 
  dplyr::filter(start < start_switch | start > end_switch)
hapsdata <- hapsdata %>% 
  dplyr::bind_rows(sim_data_bb2$haplotypes %>% 
                     dplyr::filter(start > start_switch & start <= end_switch))
hapsdata <- hapsdata %>% 
  dplyr::bind_rows(sim_data_bb4$haplotypes %>% 
                     dplyr::filter(start < start_switch | start > end_switch))
hapsdata <- hapsdata %>% 
  dplyr::bind_rows(sim_data_bb3$haplotypes %>% 
                     dplyr::filter(start > start_switch & start <= end_switch))

#all the above creates a copy number structure where there is a small segment that is LOH that is unique to a clone,
#this causes phasing issues because phasing is done relative to the clone with the largest amount of LOH per chromosome.
#to see what this looks like do plotHeatmap

results <- callHaplotypeSpecificCN(cndata %>% dplyr::filter(chr == "1"), 
                                   hapsdata %>% dplyr::filter(chr == "1"), 
                                   likelihood = "betabinomial")
results_fix <- callHaplotypeSpecificCN(cndata %>% dplyr::filter(chr == "1"), 
                                   hapsdata %>% dplyr::filter(chr == "1"), 
                                   likelihood = "betabinomial",
                                   global_phasing_for_balanced = TRUE)

nsegs_per_cell <- create_segments(results$data, field = "state_phase") %>% 
  dplyr::group_by(cell_id) %>% 
  dplyr::summarize(n = dplyr::n())

nsegs_per_cell_fix <- create_segments(results_fix$data, field = "state_phase") %>% 
  dplyr::group_by(cell_id) %>% 
  dplyr::summarize(n = dplyr::n())

test_that("Test that using all cells for phasing regions that are diploid in phasing clusters improves inference", {
  expect_true(all(nsegs_per_cell_fix$n == 3))
  expect_gt(nsegs_per_cell$n %>% mean, nsegs_per_cell_fix$n %>% mean)
})

