# Simulated joint circular-linear regression with a cyclic RW2 effect

if (!requireNamespace("INLA", quietly = TRUE)) {
  stop("This example requires the 'INLA' package.")
}

suppressPackageStartupMessages(library(INLA))
suppressPackageStartupMessages(library(INLAcircular))

set.seed(20260803)

n <- 300
n_cycle <- 20
a0 <- 1.0
a1 <- 1.5
kappa <- 300
b0 <- -0.5
b1 <- 1.2
tau <- 100
alpha <- 2
beta <- 3

z <- stats::rnorm(n, mean = 2, sd = 4)

Q <- INLA:::inla.rw2(
  n = n_cycle,
  sparse = TRUE,
  scale.model = TRUE,
  cyclic = TRUE
)
diag(Q) <- diag(Q) + 1e-5
w_base <- as.numeric(INLA::inla.qsample(
  1,
  Q,
  constr = list(A = matrix(1, 1, n_cycle), e = 0)
))

w_index <- rep(seq_len(n_cycle), length.out = n)
w <- w_base[w_index]

eta_x <- a0 + a1 * w + alpha * z
x <- rlavm(n, eta = eta_x, kappa = kappa)
y <- stats::rnorm(
  n,
  mean = b0 + b1 * eta_x + beta * z,
  sd = 1 / sqrt(tau)
)

cyclic_data <- data.frame(x = x, y = y, z = z)

cyclic_model <- list(
  likelihood(
    x ~ intercept(name = "a0", mean = 0, sd = 10) +
      covariate(z, name = "alpha", mean = 0, sd = 10) +
      f(
        w,
        model = "rw2",
        cyclic = TRUE,
        constr = TRUE,
        scale.model = TRUE,
        hyper = list(
          prec = list(
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
    y ~ intercept(name = "b0", mean = 0, sd = 10) +
      covariate(z, name = "beta", mean = 0, sd = 10) +
      covariate(
        x,
        name = "b1",
        mean = 0,
        sd = 1,
        initial = 0,
        fixed = FALSE
      ),
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
  )
)

cyclic_index <- list(
  index(
    var = "w",
    data.id = w_index,
    process.id = seq_len(n_cycle)
  )
)

cyclic_fit <- inlacc(
  model = cyclic_model,
  data = cyclic_data,
  latent.index = cyclic_index,
  metrics = TRUE,
  control.predictor = list(compute = TRUE)
)

summary(cyclic_fit, decimal = 3L)
