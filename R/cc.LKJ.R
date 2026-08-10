.LKJ_cgeneric_model <- function(n, eta, prior.u, prior.alpha) {
  dependencies <- c("graphpcor", "INLAtools")
  available <- vapply(
    dependencies,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
  missing <- dependencies[!available]
  if (length(missing)) {
    stop(
      sprintf(
        "Using model = 'iidkd_LKJ' requires %s.",
        paste(sprintf("the '%s' package", missing), collapse = " and ")
      ),
      call. = FALSE
    )
  }
  lkj_builder <- get("cgeneric_LKJ", envir = asNamespace("graphpcor"),
                     inherits = FALSE)
  INLAtools::cgeneric(
    model = lkj_builder,
    n = n,
    eta = eta,
    sigma.prior.reference = prior.u,
    sigma.prior.probability = prior.alpha,
    useINLAprecomp = FALSE
  )
}

#' Prepare LKJ covariance components for `inlacc()`
#'
#' `LKJcc()` finds every `f(..., model = "iidkd_LKJ")` term, groups terms
#' that use the same index name, constructs one `graphpcor` LKJ model per
#' group, rewrites the formulas, and creates the component and replicate
#' indices required by INLA.
#'
#' @param formulas A list of model formulas.
#' @param LKJ.eta Positive scalar LKJ shape parameter. The default is `5`.
#' @param n.obs Number of observations in each likelihood block.
#' @param latent.index Optional `index()` specification(s). For an LKJ index,
#'   `process.id = "likelihood"` means that `data.id` supplies the replicate
#'   index; the LKJ component index is generated from likelihood order.
#'
#' @return A list containing rewritten formulas, cgeneric models, per-block
#'   indices, LKJ metadata, and any unconsumed `latent.index` specifications.
#'
#' @details In every `f(index, model = "iidkd_LKJ")` term,
#'   `pc.prior.u` defaults to `1` and `pc.prior.alpha` defaults to `0.5`.
#'   These arguments may be supplied separately in each likelihood when
#'   component-specific PC priors are required. The common `LKJ.eta` value
#'   defaults to `5`.
#' @export
LKJcc <- function(formulas, LKJ.eta = 5, n.obs = NULL, latent.index = NULL) {
  if (inherits(formulas, "formula")) formulas <- list(formulas)
  if (!is.list(formulas) ||
      length(formulas) == 0L ||
      !all(vapply(formulas, inherits, logical(1), what = "formula"))) {
    stop("'formulas' must be a formula or a non-empty list of formulas.",
         call. = FALSE)
  }
  if (!is.numeric(LKJ.eta) || length(LKJ.eta) != 1L ||
      is.na(LKJ.eta) || !is.finite(LKJ.eta) || LKJ.eta <= 0) {
    stop("'LKJ.eta' must be one positive finite number.", call. = FALSE)
  }

  is_lkj_term <- function(expr) {
    is.call(expr) &&
      identical(expr[[1L]], quote(f)) &&
      "model" %in% names(expr) &&
      is.character(expr[["model"]]) &&
      length(expr[["model"]]) == 1L &&
      identical(expr[["model"]], "iidkd_LKJ")
  }

  eval_scalar <- function(expr, env, label, allow_na = FALSE) {
    value <- tryCatch(
      eval(expr, envir = env),
      error = function(e) {
        stop(sprintf("Could not evaluate '%s': %s", label, conditionMessage(e)),
             call. = FALSE)
      }
    )
    valid <- is.numeric(value) && length(value) == 1L
    if (valid) {
      valid <- if (allow_na && is.na(value)) TRUE else is.finite(value)
    }
    if (!valid) {
      stop(sprintf("'%s' must evaluate to one numeric value.", label),
           call. = FALSE)
    }
    as.numeric(value)
  }

  registry <- list()

  scan_formula <- function(expr, likelihood_id, formula_env) {
    if (!is.call(expr)) return(invisible(NULL))

    if (is_lkj_term(expr)) {
      if (length(expr) < 2L || !is.name(expr[[2L]])) {
        stop("The first argument of an iidkd_LKJ f() term must be a name, for example f(i, ...).",
             call. = FALSE)
      }

      var_name <- as.character(expr[[2L]])
      if ("eta" %in% names(expr)) {
        stop("Specify the LKJ parameter with inlacc(LKJ.eta = ...), not inside f().",
             call. = FALSE)
      }
      if ("prior.param" %in% names(expr)) {
        stop("Use 'pc.prior.u' and 'pc.prior.alpha' in each iidkd_LKJ f() term; 'prior.param' is not supported.",
             call. = FALSE)
      }

      u_value <- if ("pc.prior.u" %in% names(expr)) {
        eval_scalar(expr[["pc.prior.u"]], formula_env,
                    sprintf("pc.prior.u for '%s'", var_name))
      } else {
        1
      }
      alpha_value <- if ("pc.prior.alpha" %in% names(expr)) {
        eval_scalar(expr[["pc.prior.alpha"]], formula_env,
                    sprintf("pc.prior.alpha for '%s'", var_name),
                    allow_na = TRUE)
      } else {
        0.5
      }

      if (!is.finite(u_value) || u_value <= 0) {
        stop(sprintf("'pc.prior.u' for iidkd_LKJ process '%s' must be positive.",
                     var_name), call. = FALSE)
      }
      if (!is.na(alpha_value) && (alpha_value < 0 || alpha_value > 1)) {
        stop(sprintf("'pc.prior.alpha' for iidkd_LKJ process '%s' must be in [0, 1] or NA.",
                     var_name), call. = FALSE)
      }

      replicate_name <- NA_character_
      if ("replicate" %in% names(expr)) {
        if (!is.name(expr[["replicate"]])) {
          stop(sprintf("The 'replicate' argument for iidkd_LKJ process '%s' must be a variable name.",
                       var_name), call. = FALSE)
        }
        replicate_name <- as.character(expr[["replicate"]])
      }

      current <- registry[[var_name]]
      if (is.null(current)) {
        current <- list(
          likelihoods = integer(),
          pc.prior.u = numeric(),
          pc.prior.alpha = numeric(),
          replicate.names = character()
        )
      }
      if (likelihood_id %in% current$likelihoods) {
        stop(sprintf(
          "Likelihood %d contains more than one iidkd_LKJ term for process '%s'.",
          likelihood_id, var_name
        ), call. = FALSE)
      }

      current$likelihoods <- c(current$likelihoods, likelihood_id)
      current$pc.prior.u <- c(current$pc.prior.u, u_value)
      current$pc.prior.alpha <- c(current$pc.prior.alpha, alpha_value)
      current$replicate.names <- c(current$replicate.names, replicate_name)
      registry[[var_name]] <<- current
      return(invisible(NULL))
    }

    for (k in seq_along(expr)) {
      scan_formula(expr[[k]], likelihood_id, formula_env)
    }
    invisible(NULL)
  }

  for (i in seq_along(formulas)) {
    scan_formula(formulas[[i]], i, environment(formulas[[i]]))
  }

  empty_blocks <- lapply(seq_along(formulas), function(i) list())
  if (length(registry) == 0L) {
    return(list(
      formulas = formulas,
      models = list(),
      block.data = empty_blocks,
      lkj_vars = character(),
      replicate_vars = list(),
      likelihoods = list(),
      processes = list(),
      data.vars = character(),
      latent.index = latent.index
    ))
  }

  if (!is.numeric(n.obs) || length(n.obs) != 1L ||
      is.na(n.obs) || n.obs < 1 || n.obs != as.integer(n.obs)) {
    stop("'n.obs' must be supplied as one positive integer when iidkd_LKJ is used.",
         call. = FALSE)
  }
  n.obs <- as.integer(n.obs)

  validate_other_f_terms <- function(expr, likelihood_id) {
    if (!is.call(expr)) return(invisible(NULL))
    if (identical(expr[[1L]], quote(f)) && length(expr) >= 2L && is.name(expr[[2L]])) {
      var_name <- as.character(expr[[2L]])
      if (var_name %in% names(registry) && !is_lkj_term(expr)) {
        stop(sprintf(
          "Index '%s' is used by both iidkd_LKJ and another f() model (likelihood %d). Use different index names.",
          var_name, likelihood_id
        ), call. = FALSE)
      }
      return(invisible(NULL))
    }
    for (k in seq_along(expr)) validate_other_f_terms(expr[[k]], likelihood_id)
    invisible(NULL)
  }
  for (i in seq_along(formulas)) validate_other_f_terms(formulas[[i]], i)

  # Resolve one replicate-column name per LKJ process. Explicit names are
  # supported for backward compatibility; otherwise a collision-resistant
  # internal name is generated.
  safe_names <- stats::setNames(
    make.unique(make.names(names(registry))),
    names(registry)
  )
  for (var_name in names(registry)) {
    explicit_names <- as.character(unique(stats::na.omit(
      registry[[var_name]]$replicate.names
    )))
    if (length(explicit_names) > 1L) {
      stop(sprintf("Inconsistent replicate variables for iidkd_LKJ process '%s'.",
                   var_name), call. = FALSE)
    }
    registry[[var_name]]$replicate.var <- if (length(explicit_names) == 1L) {
      explicit_names
    } else {
      paste0(".cc_lkj_replicate_", safe_names[[var_name]])
    }
    registry[[var_name]]$model.var <-
      paste0(".cc_lkj_model_", safe_names[[var_name]])
  }

  # Consume the special LKJ index mapping. All other index() specifications
  # are returned unchanged for inlacc()'s ordinary latent-index handling.
  index_list <- latent.index
  if (!is.null(index_list) && !is.null(index_list$var)) index_list <- list(index_list)
  if (is.null(index_list)) index_list <- list()
  if (!is.list(index_list)) {
    stop("'latent.index' must be an index() object or a list of index() objects.",
         call. = FALSE)
  }

  replicate_values <- list()
  remaining_index <- list()
  for (idx in index_list) {
    if (!is.list(idx) || is.null(idx$var)) {
      stop("Every element of 'latent.index' must be created by index().",
           call. = FALSE)
    }
    idx_var <- as.character(idx$var)[1L]
    is_lkj_mapping <- idx_var %in% names(registry) &&
      identical(idx$process.id, "likelihood")

    if (!is_lkj_mapping) {
      remaining_index[[length(remaining_index) + 1L]] <- idx
      next
    }
    if (!is.null(replicate_values[[idx_var]])) {
      stop(sprintf("More than one LKJ index() mapping was supplied for '%s'.",
                   idx_var), call. = FALSE)
    }

    values <- idx$data.id
    if (!is.numeric(values) || length(values) != n.obs ||
        anyNA(values) || any(!is.finite(values)) || any(values < 1) ||
        any(values != as.integer(values))) {
      stop(sprintf(
        "For index(var = '%s', ..., process.id = 'likelihood'), data.id must contain %d positive integers.",
        idx_var, n.obs
      ), call. = FALSE)
    }
    replicate_values[[idx_var]] <- as.integer(values)
  }
  if (length(remaining_index) == 0L) remaining_index <- NULL

  models <- list()
  block_data <- empty_blocks
  processes <- list()

  add_block_column <- function(block_id, column, value) {
    old <- block_data[[block_id]][[column]]
    if (!is.null(old) && !identical(old, value)) {
      stop(sprintf("Conflicting automatically generated LKJ data for column '%s'.",
                   column), call. = FALSE)
    }
    block_data[[block_id]][[column]] <<- value
  }

  for (var_name in names(registry)) {
    info <- registry[[var_name]]
    dimension <- length(info$likelihoods)
    if (dimension < 2L) {
      stop(sprintf(
        "iidkd_LKJ process '%s' occurs in only one likelihood; an LKJ covariance model requires at least two.",
        var_name
      ), call. = FALSE)
    }

    models[[info$model.var]] <- .LKJ_cgeneric_model(
      n = dimension,
      eta = as.numeric(LKJ.eta),
      prior.u = info$pc.prior.u,
      prior.alpha = info$pc.prior.alpha
    )

    reps <- replicate_values[[var_name]]
    if (is.null(reps)) reps <- seq_len(n.obs)
    for (component in seq_along(info$likelihoods)) {
      likelihood_id <- info$likelihoods[component]
      add_block_column(likelihood_id, var_name,
                       rep.int(as.integer(component), n.obs))
      add_block_column(likelihood_id, info$replicate.var, reps)
    }

    processes[[var_name]] <- list(
      dimension = dimension,
      likelihoods = info$likelihoods,
      eta = as.numeric(LKJ.eta),
      pc.prior.u = info$pc.prior.u,
      pc.prior.alpha = info$pc.prior.alpha,
      model.var = info$model.var,
      replicate.var = info$replicate.var
    )
  }

  rewrite_formula <- function(expr) {
    if (!is.call(expr)) return(expr)
    if (is_lkj_term(expr)) {
      var_name <- as.character(expr[[2L]])
      info <- registry[[var_name]]
      expr[["model"]] <- as.name(info$model.var)
      expr[["pc.prior.u"]] <- NULL
      expr[["pc.prior.alpha"]] <- NULL
      expr[["replicate"]] <- as.name(info$replicate.var)
      return(expr)
    }
    for (k in seq_along(expr)) expr[[k]] <- rewrite_formula(expr[[k]])
    expr
  }
  formulas <- lapply(formulas, rewrite_formula)

  replicate_vars <- lapply(processes, `[[`, "replicate.var")
  likelihoods <- lapply(processes, `[[`, "likelihoods")
  data_vars <- unique(c(names(processes), unname(unlist(replicate_vars))))

  list(
    formulas = formulas,
    models = models,
    block.data = block_data,
    lkj_vars = names(processes),
    replicate_vars = replicate_vars,
    likelihoods = likelihoods,
    processes = processes,
    data.vars = data_vars,
    latent.index = remaining_index
  )
}
