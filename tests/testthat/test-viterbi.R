# Up to signals 0.16.0 the backtrace in both viterbi() and viterbiR() was not a
# Viterbi backtrace: it seeded the final position with the predecessor of the
# best final state, and then took the per-column argmax of T1 at each step
# rather than following the stored backpointers. Because it ignored the
# backpointers it also largely ignored the transition penalty, so it emitted
# isolated single-bin state changes that are absent from the MAP path.

path_loglik <- function(B, emission, transition) {
  sum(emission[cbind(seq_along(B), B + 1)]) + log(1 / nrow(transition)) +
    sum(transition[cbind(B[-length(B)] + 1, B[-1] + 1)])
}

# exhaustive search over every state path - the ground truth for small problems
best_path_bruteforce <- function(emission, transition) {
  nS <- nrow(transition); nO <- nrow(emission)
  grid <- as.matrix(expand.grid(rep(list(seq_len(nS) - 1L), nO)))
  ll <- apply(grid, 1, path_loglik, emission = emission, transition = transition)
  as.integer(grid[which.max(ll), ])
}

flat_transition <- function(nS, stp) {
  tr <- matrix(log((1 - stp) / (nS - 1)), nS, nS)
  diag(tr) <- log(stp)
  tr
}

test_that("viterbi returns the maximum a posteriori path", {
  set.seed(1)
  for (rep in 1:25) {
    nS <- 3; nO <- 6
    emission <- matrix(log(runif(nO * nS)), nO, nS)
    tr <- flat_transition(nS, runif(1, 0.6, 0.98))
    truth <- best_path_bruteforce(emission, tr)
    expect_equal(as.integer(viterbi(emission, tr, seq_len(nO))), truth)
    expect_equal(as.integer(viterbiR(emission, tr, seq_len(nO))), truth)
  }
})

test_that("the C++ and R implementations agree, in both modes", {
  set.seed(2)
  for (rep in 1:25) {
    nS <- 5; nO <- 20
    emission <- matrix(log(runif(nO * nS)), nO, nS)
    tr <- flat_transition(nS, runif(1, 0.7, 0.99))
    obs <- seq_len(nO)
    expect_equal(as.integer(viterbi(emission, tr, obs)),
                 as.integer(viterbiR(emission, tr, obs)))
    expect_equal(as.integer(viterbi(emission, tr, obs, legacy = TRUE)),
                 as.integer(viterbiR(emission, tr, obs, legacy = TRUE)))
  }
})

test_that("the legacy decode is suboptimal, and the fix is what removes that", {
  set.seed(3)
  worse <- 0
  for (rep in 1:50) {
    nS <- 4; nO <- 12
    emission <- matrix(log(runif(nO * nS)), nO, nS)
    tr <- flat_transition(nS, 0.9)
    obs <- seq_len(nO)
    fixed  <- as.integer(viterbi(emission, tr, obs))
    legacy <- as.integer(viterbi(emission, tr, obs, legacy = TRUE))
    # Viterbi is optimal, so the legacy path can never score higher
    expect_lte(path_loglik(legacy, emission, tr),
               path_loglik(fixed, emission, tr) + 1e-9)
    if (path_loglik(legacy, emission, tr) < path_loglik(fixed, emission, tr) - 1e-9) {
      worse <- worse + 1
    }
  }
  # the bug is not a rare edge case
  expect_gt(worse, 25)
})

test_that("the fix removes spurious single-bin state changes", {
  # A piecewise-constant truth with weak emissions - signals runs at 0.01-0.1X,
  # so the low-signal regime is the relevant one. The MAP path should recover
  # the two real breakpoints; extra segments in the legacy decode are artefacts.
  # Asserted over replicates rather than one seed: the size of the gap depends
  # on signal strength, but the direction never does.
  nseg <- function(x) sum(x != c(x[1], x[-length(x)]))
  simulate <- function(peak, nS = 4L, nO = 200L) {
    truth <- rep(c(0L, 2L, 1L), times = c(80, 60, 60))
    t(vapply(seq_len(nO), function(i) {
      l <- log(rep((1 - peak) / (nS - 1), nS)); l[truth[i] + 1] <- log(peak)
      l + log(runif(nS, 0.5, 1.5))
    }, numeric(nS)))
  }
  set.seed(4)
  more <- fewer <- 0
  for (rep in 1:40) {
    emission <- simulate(peak = 0.25)
    tr <- flat_transition(4L, 0.9)
    obs <- seq_len(nrow(emission))
    f <- nseg(as.integer(viterbi(emission, tr, obs)))
    l <- nseg(as.integer(viterbi(emission, tr, obs, legacy = TRUE)))
    if (l > f) more <- more + 1
    if (l < f) fewer <- fewer + 1
  }
  # the legacy decode can only ever add segments, never remove them
  expect_equal(fewer, 0)
  expect_gt(more, 30)

  # with a strong signal the MAP path is exactly the two real breakpoints
  set.seed(5)
  emission <- simulate(peak = 0.55)
  tr <- flat_transition(4L, 0.98)
  expect_equal(nseg(as.integer(viterbi(emission, tr, seq_len(nrow(emission))))), 2L)
})
