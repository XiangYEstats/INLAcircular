.INLAcircular_call_argument <- function(call, name, envir, default = NULL) {
  if (name %in% names(call)) {
    eval(call[[name]], envir = envir)
  } else {
    default
  }
}

.INLAcircular_control_family_items <- function(control.family, n.family) {
  if (is.null(control.family) || length(control.family) == 0L) {
    return(rep(list(list()), n.family))
  }
  if (!is.list(control.family)) {
    stop("'control.family' must be a list.", call. = FALSE)
  }

  if (n.family == 1L) {
    return(list(.normalize_lavm_family_setting(control.family)))
  }
  if (length(control.family) != n.family) {
    stop(
      paste0(
        "With multiple likelihoods, 'control.family' must contain one ",
        "list for each element of 'family'."
      ),
      call. = FALSE
    )
  }
  if (!all(vapply(control.family, is.list, logical(1L)))) {
    stop("Every element of 'control.family' must be a list.", call. = FALSE)
  }

  control.family
}

.INLAcircular_response_components <- function(response, n.family) {
  if (n.family == 1L) {
    return(list(response))
  }

  if (is.data.frame(response) || is.matrix(response)) {
    if (ncol(response) != n.family) {
      stop(
        "The response must have one column for each likelihood family.",
        call. = FALSE
      )
    }
    return(lapply(seq_len(n.family), function(i) response[, i]))
  }
  if (is.list(response) && length(response) == n.family) {
    return(response)
  }

  stop(
    paste0(
      "For multiple likelihoods, the formula response must be a matrix, ",
      "data frame, or list with one component per family."
    ),
    call. = FALSE
  )
}

.INLAcircular_data_list <- function(data) {
  if (is.null(data)) {
    return(list())
  }
  if (is.data.frame(data)) {
    return(as.list(data))
  }
  if (is.list(data)) {
    return(data)
  }
  if (is.environment(data)) {
    return(as.list(data, all.names = TRUE))
  }

  stop("'data' must be a data frame, list, or environment.", call. = FALSE)
}

.INLAcircular_evaluate_response <- function(formula, data, envir) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("A two-sided formula is required for family = 'lavm'.",
         call. = FALSE)
  }

  lhs <- formula[[2L]]
  if (is.null(data)) {
    return(eval(lhs, envir = envir))
  }
  if (is.environment(data)) {
    return(eval(lhs, envir = data))
  }
  eval(lhs, envir = data, enclos = envir)
}

.INLAcircular_safe_control_inla <- function(control.inla) {
  if (is.null(control.inla)) {
    control.inla <- list()
  }
  if (!is.list(control.inla)) {
    stop("'control.inla' must be a list.", call. = FALSE)
  }

  # These are the same stable defaults used by inlacc(). User values always
  # take precedence.
  if (is.null(control.inla$cmin)) {
    control.inla$cmin <- 0
  }
  if (is.null(control.inla$compute.initial.values)) {
    control.inla$compute.initial.values <- TRUE
  }
  control.inla
}

#' Run INLA with the INLAcircular LAvM likelihood
#'
#' `INLAcircular::inla()` is a transparent compatibility layer over
#' `INLA::inla()`. Calls without `family = "lavm"` are forwarded unchanged.
#' For LAvM calls, it installs the package's C `cloglike`, converts the response
#' to an INLA matrix-valued response, and evaluates the `pc.vm0` or `pc.vminf`
#' concentration prior in the same C module used by [inlacc()].
#'
#' @param ... Arguments accepted by `INLA::inla()`. For `family = "lavm"`,
#'   use `control.family = list(link = ..., hyper = list(kappa = ...))`.
#'   The single-family nested form `list(list(...))` is also accepted.
#'
#' @details Loading `INLAcircular` after `INLA` makes this function available
#'   as `inla()`. The underlying likelihood sent to INLA is `"cloglike"`, which
#'   is necessary because external packages cannot add a new compiled family or
#'   prior keyword to an already-built INLA executable. All non-LAvM arguments
#'   retain INLA's own matching and evaluation rules.
#'
#'   Within `hyper$kappa`, `initial` is on the internal `log(kappa)` scale.
#'   The defaults are `initial = 6`, `prior = "pc.vminf"`,
#'   `param = c(0.5, 0.5)`, and `fixed = FALSE`.
#'
#' @return The fitted object returned by `INLA::inla()`.
#' @export
#' @examples
#' \dontrun{
#' library(INLA)
#' library(INLAcircular)
#'
#' set.seed(1)
#' n <- 100
#' z <- rnorm(n)
#' y <- rlavm(n, eta = 1 + 1.5 * z, kappa = 100)
#'
#' fit <- inla(
#'   y ~ z,
#'   family = "lavm",
#'   data = data.frame(y = y, z = z),
#'   control.family = list(
#'     link = "inverse.tangent",
#'     hyper = list(
#'       kappa = list(
#'         initial = 6,
#'         prior = "pc.vminf",
#'         param = c(0.5, 0.5),
#'         fixed = FALSE
#'       )
#'     )
#'   )
#' )
#' }
inla <- function(...) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "Package 'INLA' is required. Install the INLA testing version first.",
      call. = FALSE
    )
  }

  caller <- parent.frame()
  original.call <- match.call(expand.dots = TRUE)
  matched.call <- match.call(
    definition = INLA::inla,
    call = original.call,
    expand.dots = TRUE
  )
  requested.family <- .INLAcircular_call_argument(
    matched.call,
    "family",
    caller,
    default = "gaussian"
  )
  requested.family <- as.character(requested.family)

  if (!any(requested.family == "lavm")) {
    original.call[[1L]] <- quote(INLA::inla)
    return(eval(original.call, envir = caller))
  }

  formula <- .INLAcircular_call_argument(
    matched.call,
    "formula",
    caller,
    default = NULL
  )
  data <- .INLAcircular_call_argument(
    matched.call,
    "data",
    caller,
    default = NULL
  )
  control.family <- .INLAcircular_call_argument(
    matched.call,
    "control.family",
    caller,
    default = list()
  )
  control.inla <- .INLAcircular_call_argument(
    matched.call,
    "control.inla",
    caller,
    default = list()
  )

  formula.environment <- environment(formula)
  if (is.null(formula.environment)) {
    formula.environment <- caller
  }

  n.family <- length(requested.family)
  control.items <- .INLAcircular_control_family_items(
    control.family,
    n.family
  )
  response <- .INLAcircular_evaluate_response(
    formula,
    data,
    formula.environment
  )
  response.items <- .INLAcircular_response_components(response, n.family)

  internal.family <- requested.family
  internal.control <- control.items
  definitions <- vector("list", n.family)
  for (i in seq_len(n.family)) {
    if (requested.family[i] == "lavm") {
      definitions[[i]] <- .make_lavm_inla_family(control.items[[i]])
      internal.family[i] <- definitions[[i]]$family
      internal.control[[i]] <- definitions[[i]]$control
      if (!inherits(response.items[[i]], "inla.mdata")) {
        response.items[[i]] <- INLA::inla.mdata(response.items[[i]])
      }
    }
  }

  data.internal <- .INLAcircular_data_list(data)
  response.name <- ".INLAcircular_lavm_response"
  while (response.name %in% names(data.internal)) {
    response.name <- paste0(response.name, "_")
  }
  data.internal[[response.name]] <- if (n.family == 1L) {
    response.items[[1L]]
  } else {
    response.items
  }

  formula.internal <- formula
  formula.internal[[2L]] <- as.name(response.name)
  environment(formula.internal) <- formula.environment

  control.internal <- if (n.family == 1L) {
    internal.control[[1L]]
  } else {
    internal.control
  }
  control.inla.internal <- .INLAcircular_safe_control_inla(control.inla)

  evaluation.environment <- new.env(parent = caller)
  evaluation.environment$.INLAcircular_formula <- formula.internal
  evaluation.environment$.INLAcircular_family <- internal.family
  evaluation.environment$.INLAcircular_data <- data.internal
  evaluation.environment$.INLAcircular_control_family <- control.internal
  evaluation.environment$.INLAcircular_control_inla <- control.inla.internal

  matched.call[[1L]] <- quote(INLA::inla)
  matched.call$formula <- quote(.INLAcircular_formula)
  matched.call$family <- quote(.INLAcircular_family)
  matched.call$data <- quote(.INLAcircular_data)
  matched.call$control.family <- quote(.INLAcircular_control_family)
  matched.call$control.inla <- quote(.INLAcircular_control_inla)

  result <- eval(matched.call, envir = evaluation.environment)
  result$INLAcircular <- list(
    requested.family = requested.family,
    internal.family = internal.family,
    prior = lapply(definitions, function(x) {
      if (is.null(x)) NULL else x$prior.label
    })
  )

  meta <- list(base_families = requested.family)
  result <- tryCatch(
    lavm.rename.inla.output(result, meta),
    error = function(e) result
  )
  result
}
