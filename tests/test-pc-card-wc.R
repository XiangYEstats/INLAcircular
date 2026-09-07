library(INLAcircular)

lambda <- 2
tolerance <- 1e-10

# Circular densities are normalized, periodic, vectorized, and C-backed.
angles <- c(-pi, -1, 0, 1, pi)
card_expected <- (1 + 2 * 0.25 * cos(angles - 0.2)) / (2 * pi)
wc_expected <- (1 - 0.7^2) /
  (2 * pi * (1 + 0.7^2 - 2 * 0.7 * cos(angles - 0.2)))

stopifnot(
  isTRUE(all.equal(
    dcardioid(angles, mu = 0.2, kappa = 0.25),
    card_expected,
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    dwrappedcauchy(angles, mu = 0.2, kappa = 0.7),
    wc_expected,
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    dcardioid(angles + 2 * pi, mu = 0.2, kappa = 0.25),
    dcardioid(angles, mu = 0.2, kappa = 0.25),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    dwrappedcauchy(angles + 2 * pi, mu = 0.2, kappa = 0.7),
    dwrappedcauchy(angles, mu = 0.2, kappa = 0.7),
    tolerance = 1e-14
  )),
  identical(dcard(0, 0, 0.25), dcardioid(0, 0, 0.25)),
  identical(dwc(0, 0, 0.7), dwrappedcauchy(0, 0, 0.7)),
  identical(dcardioid(pi, 0, 0.5), 0),
  identical(dwrappedcauchy(0, 0, 1), Inf),
  identical(dwrappedcauchy(pi, 0, 1), 0),
  identical(dcardioid(NA_real_, 0, 0.5), NA_real_),
  identical(dwrappedcauchy(NA_real_, 0, 1), NA_real_),
  is.nan(dcardioid(Inf, 0, 0.5)),
  is.nan(dwrappedcauchy(Inf, 0, 1)),
  is.nan(suppressWarnings(dcardioid(0, 0, -0.1))),
  is.nan(suppressWarnings(dwrappedcauchy(0, 0, 1.1))),
  isTRUE(all.equal(
    dcardioid(angles, 0.2, 0.25, log = TRUE),
    log(dcardioid(angles, 0.2, 0.25)),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    dwrappedcauchy(angles, 0.2, 0.7, log = TRUE),
    log(dwrappedcauchy(angles, 0.2, 0.7)),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    integrate(function(x) dcardioid(x, 0, 0.4), -pi, pi)$value,
    1,
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    integrate(function(x) dwrappedcauchy(x, 0, 0.8), -pi, pi)$value,
    1,
    tolerance = 1e-10
  ))
)

# Cardioid PC prior with uniform base: a normalized truncated exponential.
stopifnot(
  isTRUE(all.equal(
    dpc.card0(c(0, 0.25, 0.5), lambda),
    c(2.9862174220383164, 1.8927757220087252, 1.7803591104625365),
    tolerance = tolerance
  )),
  identical(dpc.card0(-1, lambda), 0),
  identical(dpc.card0(0.51, lambda), 0),
  identical(ppc.card0(-1, lambda), 0),
  identical(ppc.card0(0, lambda), 0),
  identical(ppc.card0(0.5, lambda), 1),
  identical(ppc.card0(Inf, lambda), 1),
  isTRUE(all.equal(
    ppc.card0(0.25, lambda),
    0.5951402006183594,
    tolerance = tolerance
  )),
  identical(qpc.card0(0, lambda), 0),
  identical(qpc.card0(1, lambda), 0.5),
  isTRUE(all.equal(
    qpc.card0(0.5, lambda),
    0.20172636659506137,
    tolerance = tolerance
  )),
  isTRUE(all.equal(
    dpc.card0(1e-300, lambda),
    dpc.card0(0, lambda),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    ppc.card0(0.25, lambda, log.p = TRUE),
    log(ppc.card0(0.25, lambda)),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    qpc.card0(log(0.5), lambda, log.p = TRUE),
    qpc.card0(0.5, lambda),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    integrate(
      function(x) dpc.card0(x, lambda),
      0,
      0.5,
      subdivisions = 1000L,
      rel.tol = 1e-10
    )$value,
    1,
    tolerance = 1e-8
  ))
)

card0_probabilities <- c(1e-8, 0.01, 0.2, 0.5, 0.9, 1 - 1e-10)
stopifnot(isTRUE(all.equal(
  ppc.card0(qpc.card0(card0_probabilities, lambda), lambda),
  card0_probabilities,
  tolerance = 1e-9
)))

# The lower-tail inversion must retain subnormal kappa values.
card0_tiny_probability <- .Machine$double.xmin
card0_tiny_quantile <- qpc.card0(card0_tiny_probability, lambda)
card0_tiny_linear_approximation <-
  card0_tiny_probability / dpc.card0(0, lambda)
stopifnot(
  abs(card0_tiny_quantile / card0_tiny_linear_approximation - 1) < 1e-12,
  abs(ppc.card0(card0_tiny_quantile, lambda) /
        card0_tiny_probability - 1) < 1e-12
)

# At the smallest positive rates, the truncated exponential has its
# lambda -> 0 limit: a distribution uniform in PC distance.
card0_tiny_rate <- 2 * .Machine$double.xmin * .Machine$double.eps
stopifnot(
  abs(dpc.card0(0, card0_tiny_rate) /
        1.8052399711038224 - 1) < 1e-14,
  abs(ppc.card0(0.25, card0_tiny_rate) /
        0.4589647584871319 - 1) < 1e-14,
  abs(qpc.card0(0.5, card0_tiny_rate) /
        0.27145390726974605 - 1) < 1e-12
)

# Cardioid PC prior with base kappa=0.5: continuous part plus atom at zero.
card_atom <- exp(-lambda * sqrt(log(2)))
stopifnot(
  isTRUE(all.equal(card_atom, 0.18916999526870127, tolerance = 1e-14)),
  isTRUE(all.equal(
    dpc.card(c(0, 0.25, 0.49), lambda),
    c(0.45443264077452757, 1.0445693016976794, 6.627846581454266),
    tolerance = tolerance
  )),
  identical(dpc.card(-1, lambda), 0),
  identical(dpc.card(0.51, lambda), 0),
  identical(dpc.card(0.5, lambda), Inf),
  identical(ppc.card(-1, lambda), 0),
  isTRUE(all.equal(ppc.card(0, lambda), card_atom, tolerance = 1e-14)),
  identical(ppc.card(0.5, lambda), 1),
  identical(ppc.card(Inf, lambda), 1),
  identical(qpc.card(card_atom / 2, lambda), 0),
  identical(qpc.card(card_atom, lambda), 0),
  identical(qpc.card(1, lambda), 0.5),
  isTRUE(all.equal(
    qpc.card(0.5, lambda),
    0.35395750682338145,
    tolerance = tolerance
  )),
  isTRUE(all.equal(
    ppc.card(0.25, lambda, log.p = TRUE),
    log(ppc.card(0.25, lambda)),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    qpc.card(log(0.5), lambda, log.p = TRUE),
    qpc.card(0.5, lambda),
    tolerance = 1e-14
  )),
  is.finite(dpc.card(0.5 - .Machine$double.eps, lambda)),
  isTRUE(all.equal(
    integrate(
      function(x) dpc.card(x, lambda),
      0,
      0.5,
      subdivisions = 1000L
    )$value,
    1 - card_atom,
    tolerance = 1e-7
  ))
)

card_probabilities <- c(0.2, 0.5, 0.9, 0.99, 1 - 1e-10)
stopifnot(isTRUE(all.equal(
  ppc.card(qpc.card(card_probabilities, lambda), lambda),
  card_probabilities,
  tolerance = 1e-9
)))

# Inversion immediately above the boundary atom must avoid subtracting two
# nearly equal distances.
card_probability_above_atom <- card_atom * (1 + .Machine$double.eps)
card_quantile_above_atom <- qpc.card(card_probability_above_atom, lambda)
card_first_order_quantile <- sqrt(log(2)) *
  (log(card_probability_above_atom) - log(card_atom)) / lambda
stopifnot(
  abs(card_quantile_above_atom / card_first_order_quantile - 1) < 1e-12,
  identical(
    ppc.card(card_quantile_above_atom, lambda),
    card_probability_above_atom
  )
)

# Wrapped Cauchy PC prior: analytic inversion and stable open endpoint.
stopifnot(
  isTRUE(all.equal(
    dpc.wc(c(0, 0.5, 0.9), lambda),
    c(2, 0.8503661068392163, 0.5585012996676915),
    tolerance = tolerance
  )),
  identical(dpc.wc(-1, lambda), 0),
  identical(dpc.wc(1.1, lambda), 0),
  identical(dpc.wc(1, lambda), Inf),
  identical(ppc.wc(-1, lambda), 0),
  identical(ppc.wc(0, lambda), 0),
  identical(ppc.wc(1, lambda), 1),
  identical(ppc.wc(Inf, lambda), 1),
  identical(qpc.wc(0, lambda), 0),
  identical(qpc.wc(1, lambda), 1),
  isTRUE(all.equal(
    ppc.wc(0.5, lambda),
    0.6579232126154987,
    tolerance = tolerance
  )),
  isTRUE(all.equal(
    qpc.wc(0.5, lambda),
    0.33642236016742926,
    tolerance = tolerance
  )),
  isTRUE(all.equal(
    dpc.wc(1e-300, lambda),
    lambda,
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    ppc.wc(0.5, lambda, log.p = TRUE),
    log(ppc.wc(0.5, lambda)),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    qpc.wc(log(0.5), lambda, log.p = TRUE),
    qpc.wc(0.5, lambda),
    tolerance = 1e-14
  )),
  qpc.wc(-1e-20, lambda, log.p = TRUE) < 1,
  is.finite(dpc.wc(1 - .Machine$double.eps, lambda))
)

wc_probabilities <- c(1e-8, 0.01, 0.2, 0.5, 0.9, 0.99)
wc_quantiles <- qpc.wc(wc_probabilities, lambda)
stopifnot(
  isTRUE(all.equal(
    ppc.wc(wc_quantiles, lambda),
    wc_probabilities,
    tolerance = 1e-9
  )),
  isTRUE(all.equal(
    integrate(
      function(x) dpc.wc(x, lambda),
      0,
      qpc.wc(0.99, lambda),
      subdivisions = 1000L
    )$value,
    0.99,
    tolerance = 1e-8
  ))
)

# Log probabilities and ordinary probabilities far below sqrt(DBL_MIN) must
# survive analytic inversion without cancellation or squaring underflow.
wc_tiny_log_probability <- log(1e-20)
wc_tiny_log_quantile <- qpc.wc(
  wc_tiny_log_probability,
  lambda,
  log.p = TRUE
)
wc_tiny_probability <- 1e-200
wc_tiny_quantile <- qpc.wc(wc_tiny_probability, lambda)
wc_tiny_rate <- 2 * .Machine$double.xmin * .Machine$double.eps
wc_half_distance <- sqrt(-log1p(-0.5^2))
stopifnot(
  abs(wc_tiny_log_quantile / 5e-21 - 1) < 1e-14,
  abs(ppc.wc(wc_tiny_log_quantile, lambda, log.p = TRUE) /
        wc_tiny_log_probability - 1) < 1e-14,
  abs(wc_tiny_quantile / 5e-201 - 1) < 1e-14,
  abs(ppc.wc(wc_tiny_quantile, lambda) / wc_tiny_probability - 1) < 1e-12,
  abs(ppc.wc(0.5, wc_tiny_rate, log.p = TRUE) -
        (log(wc_tiny_rate) + log(wc_half_distance))) < 1e-12
)

# RNG functions use R's RNG stream and inverse transforms in C.
set.seed(20260831)
draws_card0 <- rpc.card0(100, lambda)
set.seed(20260831)
expected_card0 <- qpc.card0(runif(100), lambda)

set.seed(20260831)
draws_card <- rpc.card(100, lambda)
set.seed(20260831)
expected_card <- qpc.card(runif(100), lambda)

set.seed(20260831)
draws_wc <- rpc.wc(100, lambda)
set.seed(20260831)
expected_wc <- qpc.wc(runif(100), lambda)

stopifnot(
  identical(draws_card0, expected_card0),
  identical(draws_card, expected_card),
  identical(draws_wc, expected_wc),
  all(draws_card0 >= 0 & draws_card0 <= 0.5),
  all(draws_card >= 0 & draws_card <= 0.5),
  any(draws_card == 0),
  all(draws_wc >= 0 & draws_wc < 1),
  length(rpc.card0(1:4, lambda)) == 4L,
  length(rpc.card(1:4, lambda)) == 4L,
  length(rpc.wc(1:4, lambda)) == 4L,
  length(rpc.card0(0, lambda)) == 0L,
  length(rpc.card(0, lambda)) == 0L,
  length(rpc.wc(0, lambda)) == 0L
)
