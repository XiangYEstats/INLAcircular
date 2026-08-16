library(INLAcircular)

parse_lavm <- getFromNamespace(
  ".parse_lavm_family_setting",
  "INLAcircular"
)
safe_control_inla <- getFromNamespace(
  ".INLAcircular_safe_control_inla",
  "INLAcircular"
)

defaults <- parse_lavm()
stopifnot(
  identical(defaults$link, "inverse.tangent"),
  identical(defaults$prior, "pc.vminf"),
  identical(defaults$u, 0.5),
  identical(defaults$alpha, 0.5),
  identical(defaults$initial, 6),
  identical(defaults$fixed, FALSE)
)

direct <- parse_lavm(list(
  link = "inverse.tangent",
  hyper = list(
    kappa = list(
      initial = 4.25,
      prior = "pc.vm0",
      param = c(0.25, 0.75),
      fixed = TRUE
    )
  )
))
nested <- parse_lavm(list(list(
  link = "inverse.tangent",
  hyper = list(
    kappa = list(
      initial = 4.25,
      prior = "pc.vm0",
      param = c(0.25, 0.75),
      fixed = TRUE
    )
  )
)))
stopifnot(
  identical(direct$link.code, 0L),
  identical(direct$prior.code, 1L),
  identical(direct$initial, 4.25),
  identical(direct$fixed, TRUE),
  identical(direct[c("link", "prior", "u", "alpha", "initial", "fixed")],
            nested[c("link", "prior", "u", "alpha", "initial", "fixed")])
)

expect_error <- function(expr, pattern) {
  message <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  stopifnot(is.character(message), grepl(pattern, message, fixed = TRUE))
}

expect_error(
  parse_lavm(list(hyper = list(kappa = list(log.initial = 6)))),
  "'log.initial' is no longer used"
)
expect_error(
  parse_lavm(list(hyper = list(prec = list(initial = 6)))),
  "Use 'hyper$kappa', not 'hyper$prec'"
)

safe_defaults <- safe_control_inla(list())
safe_override <- safe_control_inla(list(
  cmin = 0.1,
  compute.initial.values = FALSE,
  strategy = "gaussian"
))
stopifnot(
  identical(safe_defaults$cmin, 0),
  identical(safe_defaults$compute.initial.values, TRUE),
  identical(safe_override$cmin, 0.1),
  identical(safe_override$compute.initial.values, FALSE),
  identical(safe_override$strategy, "gaussian")
)
