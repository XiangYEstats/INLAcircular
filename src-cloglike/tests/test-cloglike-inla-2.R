library(INLA)
set.seed(20260819)
n <- 80L
z <- seq(-1.2, 1.2, length.out = n)
eta <- 0.4 + 0.9 * z
kappa <- 12
base_angle <- rnorm(n, sd = 1 / sqrt(kappa))
y <- 2 * atan(tan(base_angle / 2) + eta)
Y <- inla.mdata(y)

cloglike <- inla.cloglike.define(
    model = "INLAcirc_cloglike_lavm",
    shlib = NULL, ## not need as its built'in
    lavm.link = 0,
    lavm.prior = 0,
    lavm.u = 0.5,
    lavm.alpha = 0.5,
    lavm.initial.theta = log(kappa),
    lavm.fixed.theta = 1
)

fit <- inla(
    Y ~ z,
    family = "cloglike",
    data = list(Y = Y, z = z),
    control.family = list(cloglike = cloglike),
    control.inla = list(
        cmin = 0,
        compute.initial.values = TRUE
    ),
    verbose = FALSE
)

summary(fit)
