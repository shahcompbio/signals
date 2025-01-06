library(ggplot2)
library(dplyr)
library(data.table)

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
results_bb <- callHaplotypeSpecificCN(sim_data_bb$CNbins, sim_data_bb$haplotypes, likelihood = "betabinomial")
f <- table(results_bb$data[results_bb$data$chr == "1", ]$state_AS_phased)
f <- f / sum(f)

print(results_bb)

results_df <- orderdf(results_bb$data)
truth_df <- orderdf(sim_data_bb$ascn)
nbins <- length(truth_df$cell_id)
hm <- plotHeatmap(results_bb)

test_that("Test haplotype specific copy number inference (beta-binomial)", {
  expect_gt(sum(results_df$state_AS == truth_df$state_AS) / nbins, 0.99) # test accuraccy is > 99%
  expect_equal(round(results_bb$loherror, 2), loherror)
  expect_gt(results_bb$likelihood$taronesZ, 5)
  expect_true(f[["0|2"]] == 0.125 | f[["0|2"]] == 0.25)
  expect_true(f[["2|0"]] == 0.25 | f[["2|0"]] == 0.125)
  expect_true(typeof(hm) == "S4")
})

plot1 <- plotCNprofile(sim_data_bb$CNbins)
plot2 <- plotCNprofile(sim_data_bb$CNbins, xaxis_order = "bin", chrfilt = "1", maxCN = 20)
plot3 <- plotCNprofileBAF(sim_data_bb$ascn)
plot4 <- plotCNprofileBAF(sim_data_bb$ascn, xaxis_order = "bin")
plot5 <- plotBAFperstate(results_bb)
plot6 <- plotBBfit(results_bb)
plot7 <- plot_variance_state(results_bb)
plot8 <- plot_variance_state(results_bb, by_allele_specific_state = TRUE)
plot9 <- plotCNBAF(sim_data_bb$ascn)

test_that("Test plotting", {
  expect_true(is.ggplot(plot1))
  expect_true(is.ggplot(plot2))
  expect_true(is.ggplot(plot3))
  expect_true(is.ggplot(plot4))
  expect_true(is.ggplot(plot5))
  expect_true(is.ggplot(plot6))
  expect_true(is.ggplot(plot7))
  expect_true(is.ggplot(plot8))
  expect_true(is.ggplot(plot9))
})

# test rephasing by minimizing number of events
trueA <- data.frame(x = c(1, 1, 2, 2, 4, 4, 2, 2, 2, 2, 2, 2))
trueB <- data.frame(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

# scramble the CNAs
A <- data.frame(x = c(1, 1, 1, 1, 4, 4, 2, 1, 2, 1, 2, 1))
B <- data.frame(x = c(1, 1, 2, 2, 1, 1, 1, 2, 1, 2, 1, 2))


# test rephasing by minimizing number of events
trueA <- data.frame(x = c(1, 1, 2, 2, 4, 4, 2, 2, 2, 2, 2, 2))
trueB <- data.frame(x = c(1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3))

# scramble the CNAs
A <- data.frame(x = c(1, 1, 1, 1, 4, 3, 3, 3, 2, 2, 2, 2))
B <- data.frame(x = c(1, 1, 2, 2, 1, 4, 2, 2, 3, 3, 3, 3))

x <- getphase(A, B)

mymat <- cbind(as.vector(A$x), as.vector(B$x))
newA <- c()
for (i in 1:dim(mymat)[1]) {
  newA <- c(newA, mymat[i, x$phasebin[i] + 1])
}

test_that("Test rephasing by minimizing number of events", {
  expect_true(isTRUE(all.equal(trueA$x, newA)) | isTRUE(all.equal(trueB$x, newA)))
})

df <- sim_data_bb$CNbins %>% 
  group_by(cell_id, chr) %>% 
  summarize(x = sum(state != 2) / n()) %>% 
  arrange(desc(x)) %>% 
  group_by(chr) %>% 
  filter(row_number() < 5)
chr_cell_list <- split(df$cell_id, df$chr)

results_bb_2 <- callHaplotypeSpecificCN(sim_data_bb$CNbins, 
                                      sim_data_bb$haplotypes, 
                                      likelihood = "betabinomial",
                                      chr_cell_list = chr_cell_list)
f <- table(results_bb_2$data[results_bb_2$data$chr == "1", ]$state_AS_phased)
f <- f / sum(f)
f <- sort(f)
f <- f[names(f) != "1|1"]
results_df <- orderdf(results_bb_2$data)
truth_df <- orderdf(sim_data_bb$ascn)
nbins <- length(truth_df$cell_id)

test_that("Test haplotype specific copy number inference using input chr-cell list", {
  expect_gt(sum(results_df$state_AS == truth_df$state_AS) / nbins, 0.99) # test accuraccy is > 99%
  expect_equal(round(results_bb_2$loherror, 2), loherror)
  expect_gt(results_bb_2$likelihood$taronesZ, 5)
  expect_equal(f[[1]], 0.125)
  expect_equal(f[[2]], 0.25)
})
