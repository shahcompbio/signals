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

make_haplotypes <- function(nblocks = 6, cells = c("A", "B"), chr = "1",
                            bin = 5e5, blocks_per_bin = 3) {
  data.table::rbindlist(lapply(cells, function(cid) {
    data.table::data.table(
      cell_id = cid, chr = chr,
      hap_label = seq_len(nblocks) - 1L,
      start = (floor((seq_len(nblocks) - 1L) / blocks_per_bin) * bin) + 1,
      allele0 = 5L, allele1 = 1L)
  }))[, end := start + bin - 1][]
}

test_that("blocks with no counts in the selected cells are still phased", {
  h <- make_haplotypes(nblocks = 6, cells = c("A", "B"))
  # the selected cell has no coverage of blocks 4 and 5
  h <- h[!(cell_id == "A" & hap_label %in% c(4L, 5L))]
  ph <- signals:::phase_with_fallback(h, cells = "A")
  expect_equal(nrow(ph), 6L)
  expect_false(any(is.na(ph$phase)))
  expect_setequal(ph$hap_label, 0:5)
})

test_that("uncovered blocks take the all-cell majority", {
  h <- make_haplotypes(nblocks = 2, cells = c("A", "B"))
  h <- h[!(cell_id == "A" & hap_label == 1L)]
  # flip cell B's uncovered block so the fallback is distinguishable
  h[cell_id == "B" & hap_label == 1L, `:=`(allele0 = 1L, allele1 = 5L)]
  ph <- signals:::phase_with_fallback(h, cells = "A")
  expect_equal(ph[hap_label == 0L]$phase, "allele1")  # allele0 5 > allele1 1
  expect_equal(ph[hap_label == 1L]$phase, "allele0")  # from cell B
})

test_that("phase_haplotypes_bychr keeps every block in every branch", {
  h <- make_haplotypes(nblocks = 6, cells = c("A", "B"), chr = "1")
  h <- h[!(cell_id == "A" & hap_label %in% c(4L, 5L))]
  nblocks_in <- nrow(unique(h[, c("chr", "start", "end", "hap_label"), with = FALSE]))

  ph_plain <- phase_haplotypes_bychr(ascn = NULL, haplotypes = h,
                                     chrlist = list("1" = "A"),
                                     global_phasing_for_balanced = FALSE)
  expect_equal(nrow(ph_plain), nblocks_in)
  expect_false(any(is.na(ph_plain$phase)))

  ph_arm <- phase_haplotypes_bychr(ascn = NULL, haplotypes = h,
                                   chrlist = list("1p" = "A"), phasebyarm = TRUE)
  expect_gt(nrow(ph_arm), 0)
  expect_false(any(is.na(ph_arm$phase)))
})

test_that("a block spanning two bins gets a single phase", {
  h <- data.table::data.table(
    cell_id = "A", chr = "1", hap_label = 0L,
    start = c(1, 500001), end = c(500000, 1000000),
    # counts disagree bin by bin, so per-bin phasing would give two answers
    allele0 = c(1L, 9L), allele1 = c(8L, 1L))
  ph <- phase_haplotypes(h)
  expect_equal(nrow(ph), 2L)                    # both bins still reported
  expect_equal(length(unique(ph$phase)), 1L)    # with the same phase
  # pooled counts are allele0 = 10 > allele1 = 9, so the phase is allele1
  expect_true(all(ph$phase == "allele1"))
})

test_that("per-bin phasing really would have disagreed", {
  h <- data.table::data.table(
    cell_id = "A", chr = "1", hap_label = 0L,
    start = c(1, 500001), end = c(500000, 1000000),
    allele0 = c(1L, 9L), allele1 = c(8L, 1L))
  perbin <- h[, list(a0 = sum(allele0), a1 = sum(allele1)),
              by = c("chr", "start", "end", "hap_label")]
  perbin[, phase := ifelse(a0 < a1, "allele0", "allele1")]
  expect_equal(length(unique(perbin$phase)), 2L)
})

test_that("the global_phasing_for_balanced branch keeps one phase per block", {
  BIN <- 5e5; nbins <- 6; bpb <- 2; nblocks <- nbins * bpb
  cells <- c("pA", "pB", "o1", "o2")
  h <- data.table::rbindlist(lapply(cells, function(cid) data.table::data.table(
    cell_id = cid, chr = "1", hap_label = seq_len(nblocks) - 1L,
    start = (floor((seq_len(nblocks) - 1L) / bpb) * BIN) + 1,
    allele0 = 4L, allele1 = 1L)))
  h[, end := start + BIN - 1]
  # block 3 straddles a bin boundary, so it appears in two bins
  h <- rbind(h, data.table::copy(h[hap_label == 3L])[, `:=`(start = start + BIN,
                                                            end = end + BIN)])
  # blocks 10 and 11 have no coverage in the cells used to phase
  h <- h[!(cell_id %in% c("pA", "pB") & hap_label %in% c(10L, 11L))]

  ascn <- data.table::rbindlist(lapply(cells, function(cid) data.table::data.table(
    cell_id = cid, chr = "1", start = (seq_len(nbins) - 1) * BIN + 1,
    end = seq_len(nbins) * BIN, state = c(2, 2, 2, 3, 3, 3),
    A = c(1, 1, 1, 2, 2, 2), B = 1, copy = c(2, 2, 2, 3, 3, 3))))
  ascn[, `:=`(state_AS_phased = paste0(A, "|", B), state_min = pmin(A, B),
              state_phase = ifelse(A == B, "Balanced", "A-Gained"),
              state_AS = paste0(A, "|", B), LOH = ifelse(B == 0, "LOH", "NO"),
              BAF = B / state, alleleA = 10L, alleleB = 5L, totalcounts = 15L)]

  bins_in <- unique(h[, c("chr", "start", "end", "hap_label"), with = FALSE])
  ph <- phase_haplotypes_bychr(ascn = ascn, haplotypes = data.table::copy(h),
                               chrlist = list("1" = c("pA", "pB")),
                               global_phasing_for_balanced = TRUE,
                               chrs_for_global_phasing = "1")

  expect_equal(nrow(ph), nrow(bins_in))   # nothing dropped, nothing duplicated
  expect_equal(sum(duplicated(ph[, c("chr", "start", "end", "hap_label"), with = FALSE])), 0L)
  expect_false(any(is.na(ph$phase)))
  expect_equal(length(unique(ph[hap_label == 3L]$phase)), 1L)
  expect_true(all(c(10L, 11L) %in% ph$hap_label))
})
