suppressPackageStartupMessages(library(INLAcircular))

# A missing CDF input must not enter the native interpolation binary search.
# Keep NA and NaN distinct, and ensure they do not affect adjacent valid values.
queries <- c(NA_real_, NaN, -0.5, 0.5)
valid <- 3:4

check_missing <- function(actual, expected_valid) {
  stopifnot(
    identical(actual[1L], NA_real_),
    is.nan(actual[2L]),
    isTRUE(all.equal(actual[valid], expected_valid, tolerance = 1e-14))
  )
}

for (log_probability in c(FALSE, TRUE)) {
  for (strategy in c("circular", "linear")) {
    expected <- pvm(
      queries[valid], mu = 0, kappa = 2,
      strategy = strategy, log = log_probability
    )
    check_missing(
      pvm(queries, mu = 0, kappa = 2,
          strategy = strategy, log = log_probability),
      expected
    )
    check_missing(
      pvm(c(0, 0, queries[valid]), mu = c(NA_real_, NaN, 0, 0),
          kappa = 2, strategy = strategy, log = log_probability),
      expected
    )
  }

  expected <- plavm(queries[valid], eta = 0, kappa = 2,
                    log = log_probability)
  check_missing(
    plavm(queries, eta = 0, kappa = 2, log = log_probability),
    expected
  )
  check_missing(
    plavm(c(0, 0, queries[valid]), eta = c(NA_real_, NaN, 0, 0),
          kappa = 2, log = log_probability),
    expected
  )

  # Exercise the native boundary directly with a small known linear grid.
  grid_x <- c(-1, 0, 1)
  grid_y <- c(-2, 0, 2)
  expected <- stats::plogis(2 * queries[valid], log.p = log_probability)
  check_missing(
    .Call("INLAcirc_C_pvm", queries, rep(0, length(queries)),
          1L, log_probability, grid_x, grid_y, PACKAGE = "INLAcircular"),
    expected
  )
  check_missing(
    .Call("INLAcirc_C_plavm", queries, rep(0, length(queries)),
          log_probability, grid_x, grid_y, PACKAGE = "INLAcircular"),
    expected
  )
}

# A vector n requests length(n) draws, as documented; its values are ignored.
sample_sizes <- c(10, NA_real_, -1)
for (sampler in list(
  function(n) rvm(n, mu = 0.2, kappa = 2),
  function(n) rlavm(n, eta = 0.2, kappa = 2)
)) {
  set.seed(20261005)
  vector_result <- sampler(sample_sizes)
  set.seed(20261005)
  scalar_result <- sampler(length(sample_sizes))
  stopifnot(
    length(vector_result) == length(sample_sizes),
    identical(vector_result, scalar_result)
  )
}
