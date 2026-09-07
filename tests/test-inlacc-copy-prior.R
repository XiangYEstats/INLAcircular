suppressPackageStartupMessages(library(INLAcircular))

copy_decision <- getFromNamespace(
  ".inlacc_covariate_copy_decision",
  "INLAcircular"
)
is_response_term <- getFromNamespace(
  ".inlacc_is_response_term",
  "INLAcircular"
)

assert_error <- function(expr, pattern) {
  err <- tryCatch(force(expr), error = identity)
  stopifnot(
    inherits(err, "error"),
    grepl(pattern, conditionMessage(err), fixed = TRUE)
  )
}

stopifnot(
  is.null(formals(covariate)[["predictor"]]),
  isTRUE(copy_decision("x", NULL, "x")),
  identical(copy_decision("z", NULL, "x"), FALSE),
  identical(copy_decision("x", FALSE, "x"), FALSE),
  isTRUE(copy_decision("x", TRUE, "x")),
  isTRUE(is_response_term("x", "x")),
  identical(is_response_term("log(x)", "x"), FALSE)
)

assert_error(
  copy_decision("log(x)", TRUE, "x"),
  "untransformed response from another likelihood block"
)
assert_error(
  copy_decision("x", NA, "x"),
  "must be NULL or one non-missing logical value"
)
assert_error(
  copy_decision("x", c(TRUE, FALSE), "x"),
  "must be NULL or one non-missing logical value"
)

if (requireNamespace("INLA", quietly = TRUE)) {
  set.seed(20260901)
  n <- 12L
  x <- stats::rnorm(n)
  data <- data.frame(
    x = x,
    y_auto = 0.7 * x + stats::rnorm(n),
    y_on = -0.4 * x + stats::rnorm(n),
    y_off = 0.2 * x + stats::rnorm(n),
    y_bare = 0.5 * x + stats::rnorm(n),
    y_transform = log(abs(x) + 1) + stats::rnorm(n)
  )

  model <- list(
    likelihood(
      x ~ intercept(name = "alpha_x", mean = 0, sd = 2),
      family = "gaussian"
    ),
    likelihood(
      y_auto ~ covariate(x),
      family = "gaussian"
    ),
    likelihood(
      y_on ~ covariate(x, predictor = TRUE),
      family = "gaussian"
    ),
    likelihood(
      y_off ~ covariate(x, name = "observed_x", predictor = FALSE),
      family = "gaussian"
    ),
    likelihood(y_bare ~ x, family = "gaussian"),
    likelihood(y_transform ~ log(abs(x) + 1), family = "gaussian")
  )

  fit <- inlacc(
    model,
    data = data,
    control.fixed = list(
      mean = list(default = 0.25, Intercept_y_auto = 1),
      prec = list(default = 0.5, Intercept_y_auto = 2)
    ),
    verbose = FALSE
  )

  stopifnot(
    identical(fit$inlacc_meta$priors$default_mean, 0.25),
    identical(fit$inlacc_meta$priors$default_sd, 1 / sqrt(0.5)),
    identical(fit$inlacc_meta$priors$fixed_means$Intercept_y_auto, 1),
    identical(fit$inlacc_meta$priors$fixed_sds$Intercept_y_auto, 1 / sqrt(2)),
    identical(fit$inlacc_meta$priors$fixed_means$alpha_x, 0),
    identical(fit$inlacc_meta$priors$fixed_sds$alpha_x, 2),
    identical(fit$inlacc_meta$copied_covariates, "x"),
    identical(fit$inlacc_meta$fixef_map$observed_x, "observed_x")
  )

  collect_f_calls <- function(expr) {
    if (!is.call(expr)) {
      return(list())
    }
    found <- if (identical(expr[[1L]], quote(f))) list(expr) else list()
    children <- unlist(
      lapply(as.list(expr)[-1L], collect_f_calls),
      recursive = FALSE
    )
    c(found, children)
  }

  f_calls <- collect_f_calls(fit$inlacc_meta$formula[[3L]])
  copy_calls <- Filter(function(call) !is.null(call$copy), f_calls)
  copy_names <- vapply(copy_calls, function(call) as.character(call[[2L]]), "")
  expected_copy_names <- c(
    "beta_y_auto_x",
    "beta_y_on_x",
    "beta_y_bare_x"
  )

  stopifnot(identical(sort(copy_names), sort(expected_copy_names)))

  for (name in expected_copy_names) {
    call <- copy_calls[[match(name, copy_names)]]
    beta <- eval(call$hyper, envir = baseenv())$beta
    stopifnot(
      identical(beta$initial, 0),
      identical(beta$prior, "normal"),
      identical(beta$param, c(0, 0.001)),
      identical(beta$fixed, FALSE)
    )
  }

  transformed_copy_model <- list(
    likelihood(x ~ 1, family = "gaussian"),
    likelihood(
      y_auto ~ covariate(log(abs(x) + 1), predictor = TRUE),
      family = "gaussian"
    )
  )
  assert_error(
    inlacc(transformed_copy_model, data = data, verbose = FALSE),
    "untransformed response from another likelihood block"
  )
}
