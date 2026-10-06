# Allele-aware transition cost: the HMM state is the B allele copy number with
# A = total - B, so a flat matrix over B treats "B unchanged" as the cheap
# self-transition and gives "A unchanged" no standing at all.

# brute force: score every state path and return the best, so the decoder is
# checked against the definition of the MAP path rather than against itself
brute_force_pd <- function(emission, transition, tidx) {
  K <- dim(transition)[1]; T <- nrow(emission)
  paths <- as.matrix(expand.grid(rep(list(seq_len(K)), T)))
  score <- apply(paths, 1, function(p) {
    s <- -log(K) + emission[1, p[1]]
    if (T > 1) for (t in 2:T) {
      s <- s + transition[p[t - 1], p[t], tidx[t - 1]] + emission[t, p[t]]
    }
    s
  })
  list(path = paths[which.max(score), ] - 1L, score = max(score))
}

path_score_pd <- function(path, emission, transition, tidx) {
  K <- dim(transition)[1]
  s <- -log(K) + emission[1, path[1] + 1]
  if (length(path) > 1) for (t in 2:length(path)) {
    s <- s + transition[path[t - 1] + 1, path[t] + 1, tidx[t - 1]] + emission[t, path[t] + 1]
  }
  s
}

test_that("viterbi_pd returns the maximum a posteriori path", {
  set.seed(42)
  for (rep in 1:25) {
    K <- sample(2:4, 1); T <- sample(2:5, 1); M <- sample(1:3, 1)
    emission <- matrix(log(runif(T * K)), nrow = T, ncol = K)
    transition <- array(NA_real_, dim = c(K, K, M))
    for (m in seq_len(M)) {
      p <- matrix(runif(K * K), K, K)
      transition[, , m] <- log(p / rowSums(p))
    }
    tidx <- sample(seq_len(M), T - 1, replace = TRUE)
    bf <- brute_force_pd(emission, transition, tidx)
    got <- viterbi_pd(emission, transition, tidx, seq_len(T))
    # compare scores, not paths: ties are legitimate
    expect_equal(path_score_pd(got, emission, transition, tidx), bf$score,
                 tolerance = 1e-9)
  }
})

test_that("viterbiR_pd agrees with the compiled version", {
  set.seed(7)
  for (rep in 1:15) {
    K <- sample(2:5, 1); T <- sample(2:8, 1); M <- sample(1:3, 1)
    emission <- matrix(log(runif(T * K)), nrow = T, ncol = K)
    transition <- array(NA_real_, dim = c(K, K, M))
    for (m in seq_len(M)) {
      p <- matrix(runif(K * K), K, K); transition[, , m] <- log(p / rowSums(p))
    }
    tidx <- sample(seq_len(M), T - 1, replace = TRUE)
    expect_equal(
      path_score_pd(viterbiR_pd(emission, transition, tidx, seq_len(T)),
                    emission, transition, tidx),
      path_score_pd(viterbi_pd(emission, transition, tidx, seq_len(T)),
                    emission, transition, tidx), tolerance = 1e-9)
  }
})

test_that("a single slice reproduces the fixed-matrix decoder exactly", {
  set.seed(1)
  K <- 4; T <- 12
  emission <- matrix(log(runif(T * K)), nrow = T, ncol = K)
  p <- matrix(runif(K * K), K, K); tr <- log(p / rowSums(p))
  arr <- array(tr, dim = c(K, K, 1))
  expect_equal(viterbi_pd(emission, arr, rep(1L, T - 1), seq_len(T)),
               viterbi(emission, tr, seq_len(T)))
})

test_that("on constant total copy number the allele-aware matrix is the old one", {
  # dA = -dB when total CN does not change, so the cost is 0 or 2 and the matrix
  # must be identical to the flat one - this is what keeps the change confined
  # to copy number breakpoints
  minor_cn <- 0:4; K <- length(minor_cn); p <- 0.999
  tr <- signals:::allele_aware_transitions(rep(3L, 10), minor_cn, p)
  expect_equal(dim(tr$transition), c(K, K, 1))
  flat <- matrix((1 - p) / (K - 1), K, K); diag(flat) <- p
  expect_equal(exp(tr$transition[, , 1]), flat, tolerance = 1e-12)
})

test_that("across a copy number change, one-allele moves beat two-allele moves", {
  minor_cn <- 0:4; p <- 0.999
  # total CN 2 -> 3: from B = 1 (so 1|1)
  tr <- signals:::allele_aware_transitions(c(2L, 3L), minor_cn, p)
  m <- tr$transition[, , tr$tidx[1]]
  from <- which(minor_cn == 1)
  to_A_gain <- which(minor_cn == 1)   # 1|1 -> 2|1 : B same, A changes  (cost 1)
  to_B_gain <- which(minor_cn == 2)   # 1|1 -> 1|2 : A same, B changes  (cost 1)
  to_both   <- which(minor_cn == 3)   # 1|1 -> 0|3 : both change        (cost 2)
  # the two one-allele moves are now equally likely - the asymmetry is gone
  expect_equal(m[from, to_A_gain], m[from, to_B_gain], tolerance = 1e-12)
  expect_gt(m[from, to_A_gain], m[from, to_both])

  # and under the old flat matrix they were not equal: B-unchanged was the
  # self-transition and B-gain was penalised
  K <- length(minor_cn)
  flat <- matrix((1 - p) / (K - 1), K, K); diag(flat) <- p
  expect_gt(flat[from, to_A_gain], flat[from, to_B_gain])
})

test_that("2|1 -> 0|1 is preferred over 2|1 -> 1|0", {
  minor_cn <- 0:4; p <- 0.999
  # total CN 3 -> 1
  tr <- signals:::allele_aware_transitions(c(3L, 1L), minor_cn, p)
  m <- tr$transition[, , tr$tidx[1]]
  from <- which(minor_cn == 1)      # 2|1
  keepB <- which(minor_cn == 1)     # 0|1 : B unchanged, A 2->0   (cost 1)
  both  <- which(minor_cn == 0)     # 1|0 : A 2->1 and B 1->0     (cost 2)
  expect_gt(m[from, keepB], m[from, both])
})

test_that("rows are normalised and selftransitionprob = 0 gives a uniform matrix", {
  minor_cn <- 0:3
  tr <- signals:::allele_aware_transitions(c(2L, 3L, 3L, 1L), minor_cn, 0.99)
  for (m in seq_len(dim(tr$transition)[3])) {
    expect_equal(rowSums(exp(tr$transition[, , m])), rep(1, length(minor_cn)),
                 tolerance = 1e-10)
  }
  u <- signals:::allele_aware_transitions(c(2L, 3L), minor_cn, 0.0)
  expect_equal(exp(u$transition[, , 1]),
               matrix(1 / length(minor_cn), length(minor_cn), length(minor_cn)),
               tolerance = 1e-12)
})

test_that("tidx maps each transition to the right slice", {
  minor_cn <- 0:3
  bs <- c(2L, 2L, 3L, 3L, 1L)            # ds = 0, 1, 0, -2
  tr <- signals:::allele_aware_transitions(bs, minor_cn, 0.99)
  expect_equal(length(tr$tidx), length(bs) - 1L)
  expect_equal(dim(tr$transition)[3], 3L)          # ds in {-2, 0, 1}
  expect_equal(tr$tidx[1], tr$tidx[3])             # both ds = 0
  expect_false(tr$tidx[1] == tr$tidx[2])
})

test_that("HaplotypeHMM accepts allele_aware and leaves constant-CN calls alone", {
  set.seed(3)
  nbin <- 40
  binstates <- rep(2L, nbin)
  n <- rep(20L, nbin)
  x <- rbinom(nbin, n, 0.5)
  a <- HaplotypeHMM(n, x, binstates, 0:2, selftransitionprob = 0.999,
                    allele_aware = FALSE)
  b <- HaplotypeHMM(n, x, binstates, 0:2, selftransitionprob = 0.999,
                    allele_aware = TRUE)
  expect_equal(a$minorcn, b$minorcn)
})
