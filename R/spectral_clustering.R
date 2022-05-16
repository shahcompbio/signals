create_adj_matrix <- function(A, B){
  v <- c(rbind(A, B)) #interweave A and B
  m <- matrix(0, nrow = length(v), ncol = length(v))
  
  for (i in 1:length(v)){
    for (j in 1:length(v)){
      decay <- sqrt((i - j)^2)  
      decay <- abs(i - j) ^ (1/10)
      m[i,j] <- 1 / ((1 / decay) * (sqrt((v[i] - v[j])^2) + 0.0001))
    }
  }
  
  chks <- split(1:length(v), ceiling(seq_along(1:length(v)) / 2))
  for (i in 1:length(chks)){
    m[chks[[i]][1], chks[[i]][1]] <- 0
    m[chks[[i]][2], chks[[i]][2]] <- 0
    m[chks[[i]][1], chks[[i]][2]] <- 0
    m[chks[[i]][2], chks[[i]][1]] <- 0
  }
  
  return(m)
}

get_fiedler_vec <- function(x){
  L = diag(colSums(x)) - x #calculate Laplacian
  #spec <- RSpectra::eigs_sym(L, 2, which = "SM", opts = list(ncv = 100))
  spec <- RSpectra::eigs_sym(L, 2, sigma = -1e-7)
  fiedler <- spec$vectors[,1]
  sign <- fiedler > 0
  
  chks <- split(sign, ceiling(seq_along(1:length(sign)) / 2))
  myswitch <- unlist(lapply(chks, function(x) which(x)[1])) == 1
  
  as.vector(myswitch)
}

#' @export
phase_haplotypes_spectral_clustering <- function(haplotypes, CNbins){
  
  #merge haplotypes CNbins
  dat <- dplyr::inner_join(haplotypes, CNbins, by = c("cell_id", "chr", "start", "end")) %>% 
    as.data.table()
  chrs <- unique(dat$chr)
  phased_haplotypes <- data.frame()
  for (mychr in chrs){
    message(paste0("Clustering chr: ", mychr))
    dattemp  <- dat %>%
      .[chr == mychr] %>% 
      .[, unb := state %% 2 != 0] %>% 
      .[, w := fifelse(unb, 100, 1)] %>% 
      .[, list(bafA = mean(w * allele0 / totalcounts), bafB = mean(w * allele1 / totalcounts)),
        by = .(chr, start, end, hap_label)] %>% 
      .[order(chr, start, hap_label)]
    adjm <- create_adj_matrix(dattemp$bafA, dattemp$bafB)
    sw <- get_fiedler_vec(adjm)
    dattemp$sw <- sw
    dattemp <- dattemp %>% 
      dplyr::mutate(phase = fifelse(sw == TRUE, "allele1", "allele0")) %>% 
      dplyr::select(chr, start, end, hap_label, phase) %>% 
      dplyr::filter(!is.na(phase))
    phased_haplotypes <- dplyr::bind_rows(phased_haplotypes, dattemp)
  }
  return(phased_haplotypes)
}


