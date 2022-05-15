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

dplyr::inner_join(true_phase, phased_haplotypes2, by = c("chr", "start", "hap_label")) %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarize(same = sum(switch.x == switch.y) / dplyr::n()) %>% 
  dplyr::filter(chr %in% c("1", "5", "3", "17", "9"))


x1 <- sim_data_bb$haplotypes %>% 
  dplyr::filter(chr == "9") %>% 
  dplyr::distinct(chr, start, hap_label, switch) %>% 
  dplyr::pull(switch)

dat <- dplyr::inner_join(sim_data_bb$haplotypes, sim_data_bb$CNbins)

dat  <- dat %>%
  dplyr::mutate(unb = state %% 2 != 0) %>% 
  dplyr::mutate(w = ifelse(unb, 10, 1)) %>% 
  dplyr::filter(chr == "9") %>% 
  dplyr::group_by(chr, start, end, hap_label) %>% 
  dplyr::summarize(bafA = mean(w * allele0 / totalcounts), bafB = mean(w * allele1 / totalcounts))

sw <- create_adj_matrix(dat$bafA, dat$bafB) %>% get_fiedler_vec(.)
all.equal(x1, sw)

sum(x1 != sw) / length(sw)


haplotypes <- format_haplotypes_dlp(haplotypes, CNbins) %>% 
  filter_haplotypes(., 0.1)
phased_haplotypes <- phase_haplotypes_spectral_clustering(haplotypes, CNbins)
hscn <- callHaplotypeSpecificCN(CNbins, 
                                haplotypes, 
                                likelihood = "auto",
                                phased_haplotypes = phased_haplotypes)
plotHeatmap(hscn, plotcol = "state_phase")

