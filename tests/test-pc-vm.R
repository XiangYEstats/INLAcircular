library(INLAcircular)

lambda <- 2
atom <- exp(-lambda)

# Uniform-base prior: kappa = 0 is the base model.
stopifnot(
  isTRUE(all.equal(
    dpc.vm0(c(0, 1, 5), lambda),
    c(1, 0.3085593699660433, 0.01245155582179276),
    tolerance = 1e-10
  )),
  identical(dpc.vm0(-1, lambda), 0),
  identical(ppc.vm0(-1, lambda), 0),
  identical(ppc.vm0(0, lambda), 0),
  identical(ppc.vm0(Inf, lambda), 1),
  isTRUE(all.equal(
    ppc.vm0(5, lambda), 0.8842287546869297,
    tolerance = 1e-10
  )),
  identical(qpc.vm0(0, lambda), 0),
  identical(qpc.vm0(1, lambda), Inf),
  isTRUE(all.equal(
    qpc.vm0(0.5, lambda), 0.7266852927292959,
    tolerance = 1e-9
  )),
  isTRUE(all.equal(
    dpc.vm0(1e-8, lambda), 0.99999999,
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    ppc.vm0(5, lambda, log.p = TRUE),
    log(ppc.vm0(5, lambda)),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    qpc.vm0(log(0.5), lambda, log.p = TRUE),
    qpc.vm0(0.5, lambda),
    tolerance = 1e-12
  ))
)

probabilities0 <- c(0.2, 0.5, 0.9, 0.99)
stopifnot(isTRUE(all.equal(
  ppc.vm0(qpc.vm0(probabilities0, lambda), lambda),
  probabilities0,
  tolerance = 1e-9
)))

# Point-mass-base prior: the CDF includes exp(-lambda) mass at kappa = 0.
stopifnot(
  isTRUE(all.equal(
    dpc.vminf(c(0, 1, 5), lambda),
    c(0.06766764161830634, 0.10753558310742174,
      0.03696349246393501),
    tolerance = 1e-10
  )),
  identical(dpc.vminf(-1, lambda), 0),
  identical(ppc.vminf(-1, lambda), 0),
  isTRUE(all.equal(ppc.vminf(0, lambda), atom, tolerance = 1e-12)),
  identical(ppc.vminf(Inf, lambda), 1),
  isTRUE(all.equal(
    ppc.vminf(5, lambda), 0.5204586927226681,
    tolerance = 1e-10
  )),
  identical(qpc.vminf(atom, lambda), 0),
  identical(qpc.vminf(atom / 2, lambda), 0),
  isTRUE(all.equal(
    qpc.vminf(0.5, lambda), 4.484930148737734,
    tolerance = 1e-9
  )),
  isTRUE(all.equal(
    ppc.vminf(5, lambda, log.p = TRUE),
    log(ppc.vminf(5, lambda)),
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    qpc.vminf(log(0.5), lambda, log.p = TRUE),
    qpc.vminf(0.5, lambda),
    tolerance = 1e-12
  ))
)

probabilities_inf <- c(0.2, 0.5, 0.9, 0.99)
stopifnot(isTRUE(all.equal(
  ppc.vminf(qpc.vminf(probabilities_inf, lambda), lambda),
  probabilities_inf,
  tolerance = 1e-9
)))

set.seed(20260813)
draws0 <- rpc.vm0(100, lambda)
draws_inf <- rpc.vminf(100, lambda)
stopifnot(
  length(draws0) == 100L,
  length(draws_inf) == 100L,
  all(is.finite(draws0)),
  all(is.finite(draws_inf)),
  all(draws0 >= 0),
  all(draws_inf >= 0)
)
