# Developer test: New York wind joint circular-Gamma regression

if (!requireNamespace("INLA", quietly = TRUE)) {
  stop("This example requires the 'INLA' package.")
}
if (!requireNamespace("circular", quietly = TRUE)) {
  stop("This example requires the 'circular' package.")
}

suppressPackageStartupMessages(library(INLA))
suppressPackageStartupMessages(library(INLAcircular))

# Use the one-year dataset already installed with INLAcircular. This works when
# the file is sourced, run with Rscript, or evaluated line by line in RStudio.
wind_data_environment <- new.env(parent = emptyenv())
utils::data(
  "wind_newyork",
  package = "INLAcircular",
  envir = wind_data_environment
)
if (!exists("wind_newyork", envir = wind_data_environment, inherits = FALSE)) {
  stop("The package dataset 'wind_newyork' could not be loaded.", call. = FALSE)
}
wind_newyork <- wind_data_environment$wind_newyork
rm(wind_data_environment)

stopifnot(
  nrow(wind_newyork) == 8759L,
  identical(sort(unique(wind_newyork$month)), 1:12),
  all(format(wind_newyork$datetime, "%Y", tz = "UTC") == "2010")
)

# Use the first quarter of the year, as in the original application.
dat <- wind_newyork[wind_newyork$month %in% 1:3, , drop = FALSE]
direction <- as.numeric(dat$HLY.WIND.VCTDIR)
speed <- dat$HLY.WIND.VCTSPD
temperature <- dat$HLY.TEMP.NORMAL
hour_index <- dat$hour + 1L
n <- nrow(dat)

# Rotate the circular response around its sample mean and wrap to (-pi, pi].
mean_direction <- as.numeric(circular::mean.circular(
  circular::circular(direction, units = "radians")
))
x <- atan2(
  sin(direction - mean_direction),
  cos(direction - mean_direction)
)

wind_data <- data.frame(
  x = x,
  speed = speed,
  temperature = temperature
)

wind_model <- list(
  likelihood(
    x ~ intercept(name = "a0", mean = 0, sd = 1) +
      covariate(temperature, name = "a2", mean = 0, sd = 1) +
      f(
        w,
        model = "ar",
        order = 2,
        constr = TRUE,
        hyper = list(
          prec = list(
            initial = 4,
            prior = "pc.prec",
            param = c(0.5, 0.5),
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
            param = c(0.5, 0.5),
            fixed = FALSE
          )
        )
      ) +
      f(
        w_hour,
        model = "rw2",
        cyclic = TRUE,
        hyper = list(
          prec = list(
            initial = 4,
            prior = "pc.prec",
            param = c(0.5, 0.5),
            fixed = FALSE
          )
        )
      ),
    family = "lavm",
    family.setting = list(
      link = "inverse.tangent",
      hyper = list(
        kappa = list(
          initial = 6,
          prior = "pc.vminf",
          param = c(0.5, 0.99),
          fixed = FALSE
        )
      )
    )
  ),
  likelihood(
    speed ~ intercept(name = "b0", mean = 0, sd = 1) +
      covariate(temperature, name = "beta1", mean = 0, sd = 1) +
      covariate(
        x,
        name = "b1",
        mean = 0,
        sd = 1,
        initial = 0,
        fixed = FALSE
      ) +
      f(
        s,
        model = "ar",
        order = 2,
        hyper = list(
          prec = list(
            initial = 4,
            prior = "pc.prec",
            param = c(0.5, 0.5),
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
            param = c(0.5, 0.5),
            fixed = FALSE
          )
        )
      ) +
      f(
        s_hour,
        model = "rw2",
        cyclic = TRUE,
        hyper = list(
          prec = list(
            initial = 4,
            prior = "pc.prec",
            param = c(0.5, 0.5),
            fixed = FALSE
          )
        )
      ),
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
)

wind_index <- list(
  index(var = "w", data.id = seq_len(n), process.id = seq_len(n)),
  index(var = "w_hour", data.id = hour_index, process.id = 1:24),
  index(var = "s", data.id = seq_len(n), process.id = seq_len(n)),
  index(var = "s_hour", data.id = hour_index, process.id = 1:24)
)

wind_fit <- inlacc(
  model = wind_model,
  data = wind_data,
  latent.index = wind_index,
  metrics = TRUE,
  control.predictor = list(compute = TRUE)
)

summary(wind_fit, decimal = 3L)

# Block-specific predictive measures. The circular response occupies the first
# n rows of the joint response and wind speed occupies the next n rows.
idx_direction <- seq_len(n)
idx_speed <- n + seq_len(n)

wind_predictive_metrics <- data.frame(
  response = c("wind direction", "wind speed"),
  rmse = c(
    sqrt(mean((x - wind_fit$summary.fitted.values$mean[idx_direction])^2)),
    sqrt(mean((speed - wind_fit$summary.fitted.values$mean[idx_speed])^2))
  ),
  lpml = c(
    sum(log(wind_fit$cpo$cpo[idx_direction]), na.rm = TRUE),
    sum(log(wind_fit$cpo$cpo[idx_speed]), na.rm = TRUE)
  )
)

wind_predictive_metrics
