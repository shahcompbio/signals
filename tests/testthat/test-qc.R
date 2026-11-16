# qc_per_cell() / filtercn() on inputs that are not whole-genome.

test_that("qc_per_cell counts segments correctly on a chromosome subset", {
  # nsegments used to be n() - 22, which is negative for every cell when fewer
  # than 22 autosomes are present, so filtercn() discarded all of them
  sim <- simulate_data_cohort(clone_num = c(15), coverage = 50,
                              clonal_events = list(list("1" = c(2, 0))), loherror = 0.02)
  hscn <- callHaplotypeSpecificCN(sim$CNbins, sim$haplotypes, progressbar = FALSE)

  sub <- hscn
  sub$data <- dplyr::filter(hscn$data, chr %in% c("1", "2", "3"))
  part <- qc_per_cell(sub)
  expect_true(all(part$nsegments >= 0))
  expect_gt(nrow(filtercn(sub)$qc_per_cell), 0)

  # on whole-genome data the value is exactly what n() - 22 gave
  full <- qc_per_cell(hscn)
  old <- create_segments(hscn$data) %>%
    dplyr::filter(!chr %in% c("X", "Y")) %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarize(nsegments_old = dplyr::n() - 22)
  m <- dplyr::inner_join(full, old, by = "cell_id")
  expect_equal(nrow(m), nrow(full))
  expect_equal(m$nsegments, m$nsegments_old)
})
