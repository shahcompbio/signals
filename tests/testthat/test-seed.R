library(dplyr)
library(data.table)

# The phasing path has three stochastic steps: the subsampling in min_cells()
# that sets the cluster size, the UMAP embedding used to choose which cells
# phase each chromosome, and the subsampling in fitBB(). Left unseeded these
# make repeated runs on identical input select different phasing cells, which
# can flip A/B assignments for a whole chromosome.

data(CNbins)
data(haplotypes)
haps <- format_haplotypes_dlp(haplotypes, CNbins)

test_that("min_cells is reproducible when seeded", {
  h <- as.data.table(haps)
  expect_identical(
    min_cells(h, mincells = 2, seed = 42)$ncells_forclustering,
    min_cells(h, mincells = 2, seed = 42)$ncells_forclustering
  )
  expect_identical(
    min_cells(h, mincells = 2, seed = 7)$prop,
    min_cells(h, mincells = 2, seed = 7)$prop
  )
})

test_that("umap_clustering is reproducible when seeded", {
  skip_if_not_installed("uwot")
  skip_if_not_installed("dbscan")
  a <- umap_clustering(CNbins, minPts = 5, seed = 42)$clustering
  b <- umap_clustering(CNbins, minPts = 5, seed = 42)$clustering
  a <- a[order(a$cell_id), ]
  b <- b[order(b$cell_id), ]
  expect_equal(a$clone_id, b$clone_id)
  expect_equal(a$umap1, b$umap1)
  expect_equal(a$umap2, b$umap2)
})

test_that("callHaplotypeSpecificCN is reproducible when seeded", {
  skip_if_not_installed("uwot")
  skip_if_not_installed("dbscan")
  ord <- function(x) x$data[order(x$data$cell_id, x$data$chr, x$data$start), ]

  r1 <- callHaplotypeSpecificCN(CNbins, haps, mincells = 2, progressbar = FALSE, seed = 42)
  r2 <- callHaplotypeSpecificCN(CNbins, haps, mincells = 2, progressbar = FALSE, seed = 42)

  # the cells chosen to phase each chromosome are the stochastic part, so pin
  # those as well as the downstream calls
  expect_equal(lapply(r1$phasing, sort), lapply(r2$phasing, sort))
  expect_equal(r1$haplotype_phasing, r2$haplotype_phasing)

  d1 <- ord(r1); d2 <- ord(r2)
  expect_equal(d1$state_AS_phased, d2$state_AS_phased)
  expect_equal(d1$A, d2$A)
  expect_equal(d1$B, d2$B)
})
