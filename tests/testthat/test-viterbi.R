# Regression tests for the Viterbi decoders (viterbiR and the C++ viterbi()).
#
# The decoders must return the single most likely state path. For small HMMs
# that path can be found exactly by exhaustive search over all state sequences,
# giving an independent ground truth to check against. This pins down the
# backtrace, which previously recomputed a per-column argmax instead of
# following the stored backpointer chain and so returned suboptimal paths.

# Exhaustive-search reference decoder.
# emission: numObs x numStates matrix of log emission probabilities (same
# orientation the package passes to viterbiR/viterbi). transition[a, b] is the
# log probability of moving from state a to state b. Initial state is uniform.
# Returns the 0-based most likely path.
brute_force_viterbi <- function(emission, transition) {
  numObs <- nrow(emission)
  numStates <- ncol(emission)
  log_init <- log(rep(1 / numStates, numStates))
  paths <- expand.grid(rep(list(seq_len(numStates)), numObs))
  best_score <- -Inf
  best_path <- NULL
  for (r in seq_len(nrow(paths))) {
    z <- as.integer(paths[r, ])
    score <- log_init[z[1]] + emission[1, z[1]]
    if (numObs > 1) {
      for (j in 2:numObs) {
        score <- score + transition[z[j - 1], z[j]] + emission[j, z[j]]
      }
    }
    if (score > best_score) {
      best_score <- score
      best_path <- z
    }
  }
  best_path - 1L
}

test_that("viterbiR matches exhaustive-search optimum on random small HMMs", {
  set.seed(101)
  for (trial in seq_len(200)) {
    numStates <- sample(2:4, 1)
    numObs <- sample(1:7, 1)
    emission <- matrix(log(runif(numObs * numStates)), numObs, numStates)
    transition <- matrix(log(runif(numStates * numStates)), numStates, numStates)
    obs <- seq_len(numObs)

    expect_equal(
      viterbiR(emission, transition, obs),
      brute_force_viterbi(emission, transition)
    )
  }
})

test_that("C++ viterbi matches exhaustive-search optimum on random small HMMs", {
  set.seed(202)
  for (trial in seq_len(200)) {
    numStates <- sample(2:4, 1)
    numObs <- sample(1:7, 1)
    emission <- matrix(log(runif(numObs * numStates)), numObs, numStates)
    transition <- matrix(log(runif(numStates * numStates)), numStates, numStates)
    obs <- seq_len(numObs)

    expect_equal(
      as.integer(viterbi(emission, transition, obs)),
      brute_force_viterbi(emission, transition)
    )
  }
})

test_that("R and C++ decoders agree with each other", {
  set.seed(303)
  for (trial in seq_len(100)) {
    numStates <- sample(2:5, 1)
    numObs <- sample(1:12, 1)
    emission <- matrix(log(runif(numObs * numStates)), numObs, numStates)
    transition <- matrix(log(runif(numStates * numStates)), numStates, numStates)
    obs <- seq_len(numObs)

    expect_equal(
      viterbiR(emission, transition, obs),
      as.integer(viterbi(emission, transition, obs))
    )
  }
})

test_that("single-bin sequences decode to the highest-emission state, not state 0", {
  # A one-bin chromosome must return argmax of the emission, exercising the
  # numObs == 1 edge case in both decoders.
  transition <- matrix(log(c(0.7, 0.3, 0.3, 0.7)), 2, 2)
  emission <- matrix(log(c(0.1, 0.9)), nrow = 1) # state 2 far more likely

  expect_equal(viterbiR(emission, transition, 1L), 1L)
  expect_equal(as.integer(viterbi(emission, transition, 1L)), 1L)
})

test_that("the decoded path is never lower probability than the optimum", {
  # Guards against any future regression that returns a valid-looking but
  # suboptimal path: score the decoded path and compare to the exhaustive best.
  path_logprob <- function(z0, emission, transition) {
    z <- z0 + 1L
    numObs <- nrow(emission)
    numStates <- ncol(emission)
    log_init <- log(rep(1 / numStates, numStates))
    score <- log_init[z[1]] + emission[1, z[1]]
    if (numObs > 1) {
      for (j in 2:numObs) {
        score <- score + transition[z[j - 1], z[j]] + emission[j, z[j]]
      }
    }
    score
  }

  set.seed(404)
  for (trial in seq_len(100)) {
    numStates <- sample(2:4, 1)
    numObs <- sample(2:7, 1)
    emission <- matrix(log(runif(numObs * numStates)), numObs, numStates)
    transition <- matrix(log(runif(numStates * numStates)), numStates, numStates)
    obs <- seq_len(numObs)

    optimum <- brute_force_viterbi(emission, transition)
    decoded <- viterbiR(emission, transition, obs)
    expect_gte(
      path_logprob(decoded, emission, transition),
      path_logprob(optimum, emission, transition) - 1e-9
    )
  }
})
