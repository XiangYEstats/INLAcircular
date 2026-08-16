if (requireNamespace("INLA", quietly = TRUE)) {
  suppressPackageStartupMessages(library(INLA))
  suppressPackageStartupMessages(library(INLAcircular))

  stopifnot(identical(inla, INLAcircular::inla))

  set.seed(20260816)
  n <- 60
  z <- rnorm(n)
  y <- rlavm(n, eta = 0.75 + 1.25 * z, kappa = 50)

  fit <- inla(
    y ~ z,
    family = "lavm",
    data = data.frame(y = y, z = z),
    control.family = list(
      list(
        link = "inverse.tangent",
        hyper = list(
          kappa = list(
            initial = log(50),
            prior = "pc.vminf",
            param = c(0.5, 0.5),
            fixed = FALSE
          )
        )
      )
    )
  )

  stopifnot(
    inherits(fit, "inla"),
    identical(fit$INLAcircular$requested.family, "lavm"),
    identical(fit$INLAcircular$internal.family, "cloglike"),
    identical(fit$INLAcircular$prior[[1L]], "pc.vminf(0.5, 0.5)"),
    all(c("(Intercept)", "z") %in% rownames(fit$summary.fixed)),
    "kappa for lavm observations" %in% rownames(fit$summary.hyperpar)
  )
}
