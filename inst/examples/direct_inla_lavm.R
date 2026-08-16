# Direct use of the INLAcircular LAvM likelihood through inla().
# Load INLA first so INLAcircular::inla masks INLA::inla on the search path.
library(INLA)
library(INLAcircular)

set.seed(20260816)
n <- 100
beta0 <- 1
beta1 <- 1.5
kappa <- 100
z <- rnorm(n)
y <- rlavm(n, eta = beta0 + beta1 * z, kappa = kappa)
dat <- data.frame(y = y, z = z)

# The LAvM-safe defaults are inverse.tangent, initial log(kappa) = 6,
# pc.vminf with param = c(0.5, 0.5), and fixed = FALSE.
fit_default <- inla(
  y ~ z,
  family = "lavm",
  data = dat
)

# The same model with the complete INLA-style control.family specification.
# For LAvM, initial is already on the internal log(kappa) scale.
fit_controlled <- inla(
  y ~ z,
  family = "lavm",
  data = dat,
  control.family = list(
    list(
      link = "inverse.tangent",
      hyper = list(
        kappa = list(
          initial = 6,
          prior = "pc.vminf",
          param = c(0.5, 0.5),
          fixed = FALSE
        )
      )
    )
  )
)

summary(fit_default)
summary(fit_controlled)
