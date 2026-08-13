# Simulated multivariate joint model with mixed response families

if (!requireNamespace("INLA", quietly = TRUE)) {
  stop("This example requires the 'INLA' package.")
}

suppressPackageStartupMessages(library(INLA))
suppressPackageStartupMessages(library(INLAcircular))

set.seed(20260803)

n <- 300
kappa1 <- 500
kappa2 <- 200
tau <- 30

# Fixed effects used for simulation.
a10 <- 0.2
a11 <- 0.5
alpha11 <- 0.6
alpha12 <- 0
alpha13 <- 0.5

a20 <- -0.1
a21 <- 0.2
alpha21 <- 0
alpha22 <- 0.5
alpha23 <- 0

b10 <- 1
b11 <- -0.3
b12 <- 5
beta11 <- 1
beta12 <- 1
beta13 <- 1.5

b20 <- 2
b21 <- 0.8
b22 <- 0
beta21 <- 0.4
beta22 <- 0
beta23 <- 0

z1 <- stats::rnorm(n)
z2 <- stats::rgamma(n, shape = 1, rate = 1)
z3 <- stats::rpois(n, lambda = 1)

Q1 <- INLA:::inla.rw2(
  n = n,
  sparse = TRUE,
  scale.model = TRUE,
  cyclic = FALSE
)
diag(Q1) <- diag(Q1) + 1e-5
w1 <- as.numeric(INLA::inla.qsample(
  1,
  Q1,
  constr = list(
    A = matrix(rbind(1, seq_len(n)), 2, n),
    e = c(0, 0)
  )
))

w2 <- as.numeric(scale(stats::arima.sim(
  n = n,
  model = list(ar = c(0.9, -0.3))
)))
s <- as.numeric(scale(stats::rnorm(n)))

eta_x1 <- a10 + a11 * w1 + alpha11 * z1 + alpha12 * z2 + alpha13 * z3
eta_x2 <- a20 + a21 * w2 + alpha21 * z1 + alpha22 * z2 + alpha23 * z3

x1 <- rlavm(n, eta = eta_x1, kappa = kappa1)
x2 <- rlavm(n, eta = eta_x2, kappa = kappa2)

y1 <- stats::rnorm(
  n,
  mean = b10 + b11 * eta_x1 + b12 * eta_x2 +
    beta11 * z1 + beta12 * z2 + beta13 * z3,
  sd = 1 / sqrt(tau)
)
y2 <- stats::rpois(
  n,
  lambda = exp(
    b20 + b21 * eta_x1 + b22 * eta_x2 + s +
      beta21 * z1 + beta22 * z2 + beta23 * z3
  )
)

multivariate_data <- data.frame(
  x1 = x1,
  x2 = x2,
  y1 = y1,
  y2 = y2,
  z1 = z1,
  z2 = z2,
  z3 = z3
)

multivariate_model <- list(
  likelihood(
    x1 ~ intercept(name = "a10", mean = 0, sd = 10) +
      covariate(z1, name = "alpha11", mean = 0, sd = 10) +
      covariate(z2, name = "alpha12", mean = 0, sd = 10) +
      covariate(z3, name = "alpha13", mean = 0, sd = 10) +
      f(
        w1,
        model = "rw2",
        constr = TRUE,
        scale.model = TRUE,
        hyper = list(
          prec = list(
            initial = 4,
            prior = "pc.prec",
            param = c(1, 0.5),
            fixed = FALSE
          )
        )
      ),
    family = "lavm",
    family.setting = list(
      link = "inverse.tangent",
      hyper = list(
        kappa = list(
          log.initial = 6,
          prior = "pc.vminf",
          param = c(0.5, 0.99),
          fixed = FALSE
        )
      )
    )
  ),
  likelihood(
    x2 ~ intercept(name = "a20", mean = 0, sd = 10) +
      covariate(z1, name = "alpha21", mean = 0, sd = 10) +
      covariate(z2, name = "alpha22", mean = 0, sd = 10) +
      covariate(z3, name = "alpha23", mean = 0, sd = 10) +
      f(
        w2,
        model = "ar",
        order = 2,
        constr = TRUE,
        hyper = list(
          prec = list(
            initial = 4,
            prior = "pc.prec",
            param = c(1, 0.5),
            fixed = FALSE
          ),
          pacf1 = list(
            initial = 0,
            prior = "pc.cor0",
            param = c(0.5, 0.5),
            fixed = FALSE
          ),
          pacf2 = list(
            initial = 0,
            prior = "pc.cor0",
            param = c(0.5, 0.4),
            fixed = FALSE
          )
        )
      ),
    family = "lavm",
    family.setting = list(
      link = "inverse.tangent",
      hyper = list(
        kappa = list(
          log.initial = 6,
          prior = "pc.vminf",
          param = c(0.5, 0.99),
          fixed = FALSE
        )
      )
    )
  ),
  likelihood(
    y1 ~ intercept(name = "b10") +
      covariate(z1, name = "beta11", mean = 0, sd = 10) +
      covariate(z2, name = "beta12", mean = 0, sd = 10) +
      covariate(z3, name = "beta13", mean = 0, sd = 10) +
      covariate(x1, name = "b11", mean = 0, sd = 1) +
      covariate(x2, name = "b12", mean = 0, sd = 1),
    family = "gaussian",
    family.setting = list(
      hyper = list(
        prec = list(
          initial = 6,
          prior = "pc.prec",
          param = c(1, 0.5),
          fixed = FALSE
        )
      )
    )
  ),
  likelihood(
    y2 ~ intercept(name = "b20") +
      covariate(z1, name = "beta21", mean = 0, sd = 10) +
      covariate(z2, name = "beta22", mean = 0, sd = 10) +
      covariate(z3, name = "beta23", mean = 0, sd = 10) +
      covariate(x1, name = "b21", mean = 0, sd = 1) +
      covariate(x2, name = "b22", mean = 0, sd = 1) +
      f(
        s,
        model = "iid",
        constr = TRUE,
        hyper = list(
          prec = list(
            initial = 4,
            prior = "pc.prec",
            param = c(1, 0.5),
            fixed = FALSE
          )
        )
      ),
    family = "poisson",
    family.setting = list()
  )
)

multivariate_fit <- inlacc(
  model = multivariate_model,
  data = multivariate_data,
  metrics = TRUE,
  control.predictor = list(compute = TRUE)
)

summary(multivariate_fit, decimal = 3L)
