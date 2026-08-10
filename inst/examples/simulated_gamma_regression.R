# Simulated Gamma regression

if (!requireNamespace("INLA", quietly = TRUE)) {
  stop("This example requires the 'INLA' package.")
}

suppressPackageStartupMessages(library(INLA))
suppressPackageStartupMessages(library(INLAcircular))

set.seed(20260803)

n <- 300
precision <- 100
beta0 <- 2
beta1 <- 3

z <- stats::rnorm(n, mean = 2, sd = 4)
eta_y <- beta0 + beta1 * z
y <- stats::rgamma(
  n,
  shape = precision,
  rate = precision / exp(eta_y)
)

gamma_data <- data.frame(y = y, z = z)

gamma_model <- likelihood(
  y ~ intercept(name = "beta0", mean = 0, sd = 10) +
    covariate(z, name = "beta1", mean = 0, sd = 10),
  family = "gamma",
  family.setting = list(
    link = "log",
    hyper = list(
      prec = list(
        initial = 6,
        prior = "loggamma",
        param = c(1, 0.01),
        fixed = FALSE
      )
    )
  )
)

gamma_fit <- inlacc(
  model = gamma_model,
  data = gamma_data,
  metrics = TRUE,
  control.predictor = list(compute = TRUE)
)

summary(gamma_fit, decimal = 3L)
