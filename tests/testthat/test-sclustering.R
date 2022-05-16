library(ggplot2)
library(dplyr)
library(data.table)

loherror <- 0.02
sim_data_bb <- simulate_data_cohort(
  clone_num = c(20, 25, 25, 1),
  clonal_events = list(
    list("1" = c(2, 0), "5" = c(3, 1)),
    list("2" = c(6, 3), "3" = c(1, 0)),
    list("17" = c(3, 1), "8" = c(6, 2)),
    list("1" = c(2, 2), "9" = c(3, 1))
  ), # opposite LOH on chr 1
  loherror = loherror,
  coverage = 100,
  rho = 0.02,
  likelihood = "betabinomial",
  nchr = 0
)

phased_haplotypes2 <- phase_haplotypes_spectral_clustering(sim_data_bb$haplotypes, sim_data_bb$CNbins) %>% 
  dplyr::mutate(switch = ifelse(phase == "allele1", TRUE, FALSE))
true_phase <- sim_data_bb$haplotypes %>% 
  dplyr::distinct(chr, start, hap_label, switch)

test_chr <- dplyr::inner_join(true_phase, phased_haplotypes2, by = c("chr", "start", "hap_label")) %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(same = abs(cor(switch.x,switch.y))) %>% #correlation should be ~-1 or 1
  dplyr::filter(chr %in% c("1", "5", "3", "17", "9")) #chromosomes with imbalances
test_chr

test_that("Test that phasing using spectral clustering is accurate", {
  expect_true(all(test_chr$same > 0.9))
})

haplotypes <- format_haplotypes_dlp(haplotypes, CNbins) %>% 
  filter_haplotypes(., 0.1)
phased_haplotypes <- phase_haplotypes_spectral_clustering(haplotypes %>% dplyr::filter(chr == 9), CNbins)
hscn <- callHaplotypeSpecificCN(CNbins, 
                                haplotypes, 
                                likelihood = "binomial",
                                phased_haplotypes = phased_haplotypes)
plotHeatmap(hscn, plotcol = "state_phase")

v <- seq(1,100,1)
m <- c()
i<-50
for (j in 1:length(v)){
  decay <- log10(abs(i - j)^1 + 1)  
  decay <- sqrt((i - j)^2)
  decay <- abs(i - j) ^ (1/10)
  m[j]<- (1 / (decay + 1))
}
m
plot(m)
