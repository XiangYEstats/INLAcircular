# Simulated joint circular-linear regression with an RW2 effect

if (!requireNamespace("INLA", quietly = TRUE)) {
  stop("This example requires the 'INLA' package.")
}

suppressPackageStartupMessages(library(INLA))
suppressPackageStartupMessages(library(INLAcircular))

set.seed(20260803)

# Simulate a circular response x and a Gaussian response y that share the
# latent predictor of x.
n <- 300
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
  n = n,
  sparse = TRUE,
  scale.model = TRUE,
  cyclic = FALSE
)
diag(Q) <- diag(Q) + 1e-5
w <- as.numeric(INLA::inla.qsample(
  1,
  Q,
  constr = list(
    A = matrix(rbind(1, seq_len(n)), 2, n),
    e = c(0, 0)
  )
))

eta_x <- a0 + a1 * w + alpha * z
x <- rlavm(n, eta = eta_x, kappa = kappa)
y <- stats::rnorm(
  n,
  mean = b0 + b1 * eta_x + beta * z,
  sd = 1 / sqrt(tau)
)

joint_data <- data.frame(x = x, y = y, z = z)

joint_model <- list(
  likelihood(
    x ~ intercept(name = "a0", mean = 0, sd = 10) +
      covariate(z, name = "alpha", mean = 0, sd = 10) +
      f(
        w,
        model = "rw2",
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
          prior = "pc.vm.inf",
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

joint_fit <- inlacc(
  model = joint_model,
  data = joint_data,
  metrics = TRUE,
  control.predictor = list(compute = TRUE)
)

summary(joint_fit, decimal = 3L)
