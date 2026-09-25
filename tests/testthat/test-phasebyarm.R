library(dplyr)
library(data.table)

# phase_haplotypes_bychr(phasebyarm = TRUE) filters `chrarm == names(chrlist)[i]`,
# so the cell list must be keyed by chromosome arm. get_cells_per_chr_local
# previously accepted phasebyarm and ignored it, keying by chromosome instead,
# which matches nothing ("4p" != "4") and silently produced an empty phasing
# table on the default cluster_per_chr = TRUE path.

# ncells_for_clustering = 1 takes the branch that gives every cell its own
# clone, so these tests exercise the keying without needing uwot/dbscan.
make_chr4 <- function(ncells = 6L) {
  cells <- paste0("c", seq_len(ncells))
  starts <- c(seq(1, by = 5e5, length.out = 20),            # p arm
              seq(60000001, by = 5e5, length.out = 20))     # q arm
  ascn <- data.table::CJ(cell_id = cells, start = starts)
  ascn[, `:=`(chr = "4", end = start + 5e5 - 1, state = 2L, LOH = "NO",
              state_BAF = 0.5,
              balance = as.integer(start < 49660117))]      # imbalanced on p only
  haps <- data.table::CJ(cell_id = cells, start = starts)
  haps[, `:=`(chr = "4", end = start + 5e5 - 1,
              hap_label = rep(seq_along(starts), ncells),
              allele1 = 7L, allele0 = 2L)]
  list(ascn = ascn, haps = haps)
}

test_that("get_cells_per_chr_local keys by arm when phasebyarm = TRUE", {
  d <- make_chr4()
  by_chr <- get_cells_per_chr_local(d$ascn, d$haps, ncells_for_clustering = 1,
                                    phasebyarm = FALSE)
  by_arm <- get_cells_per_chr_local(d$ascn, d$haps, ncells_for_clustering = 1,
                                    phasebyarm = TRUE)
  expect_equal(sort(names(by_chr)), "4")
  expect_equal(sort(names(by_arm)), c("4p", "4q"))
})

test_that("the arm-keyed list actually phases blocks", {
  d <- make_chr4()
  for (pba in c(FALSE, TRUE)) {
    cl <- get_cells_per_chr_local(d$ascn, d$haps, ncells_for_clustering = 1,
                                  phasebyarm = pba)
    ph <- phase_haplotypes_bychr(ascn = d$ascn, haplotypes = d$haps, chrlist = cl,
                                 phasebyarm = pba,
                                 global_phasing_for_balanced = FALSE)
    # every block must receive a phase, under either unit
    expect_equal(nrow(ph), data.table::uniqueN(d$haps$hap_label))
    expect_true(all(ph$phase %in% c("allele0", "allele1")))
  }
})

test_that("a chromosome-keyed list phases nothing when phasing by arm", {
  # the regression this fixes: keys and filter must agree
  d <- make_chr4()
  mismatched <- get_cells_per_chr_local(d$ascn, d$haps, ncells_for_clustering = 1,
                                        phasebyarm = FALSE)   # keyed "4"
  ph <- phase_haplotypes_bychr(ascn = d$ascn, haplotypes = d$haps,
                               chrlist = mismatched, phasebyarm = TRUE,
                               global_phasing_for_balanced = FALSE)
  expect_equal(nrow(ph), 0)
})

test_that("proportion_imbalance forwards phasebyarm to the selection step", {
  # it previously accepted the argument and dropped it, so neither the local nor
  # the global path ever saw it
  expect_true("phasebyarm" %in% names(formals(get_cells_per_chr_local)))
  expect_true("phasebyarm" %in% names(formals(get_cells_per_chr_global)))
  body_txt <- paste(deparse(body(proportion_imbalance)), collapse = " ")
  expect_match(body_txt, "phasebyarm = phasebyarm")
})
