.canonical_lavm_pc_prior <- function(prior) {
  prior <- as.character(prior)
  if (length(prior) != 1L || is.na(prior) ||
      !(prior %in% c("pc.vm0", "pc.vminf"))) {
    stop(
      "The LAvM concentration prior must be 'pc.vm0' or 'pc.vminf'.",
      call. = FALSE
    )
  }
  prior
}

.normalize_lavm_family_setting <- function(family.setting = NULL) {
  if (is.null(family.setting)) {
    return(list())
  }
  if (!is.list(family.setting)) {
    stop("The LAvM family setting must be a list.", call. = FALSE)
  }

  # INLA accepts a single-family control either directly or as
  # list(list(...)).  Accept both forms so the direct inla() interface follows
  # the usual single- and multiple-likelihood conventions.
  outer_names <- names(family.setting)
  outer_is_unnamed <- is.null(outer_names) ||
    all(is.na(outer_names) | outer_names == "")
  if (length(family.setting) == 1L && outer_is_unnamed &&
      is.list(family.setting[[1L]])) {
    family.setting <- family.setting[[1L]]
  }

  family.setting
}

.parse_lavm_family_setting <- function(family.setting = NULL) {
  family.setting <- .normalize_lavm_family_setting(family.setting)

  link <- if (is.null(family.setting$link)) {
    "inverse.tangent"
  } else {
    as.character(family.setting$link)
  }
  if (length(link) != 1L || is.na(link)) {
    stop("The LAvM link must be one character string.", call. = FALSE)
  }
  link.code <- switch(
    link,
    "inverse.tangent" = 0L,
    "scaled.logit" = 1L,
    "scaled.probit" = 2L,
    NULL
  )
  if (is.null(link.code)) {
    stop(
      paste0(
        "Invalid LAvM link '", link, "'. Valid choices are ",
        "'inverse.tangent', 'scaled.logit', and 'scaled.probit'."
      ),
      call. = FALSE
    )
  }

  hyper <- family.setting$hyper
  if (is.null(hyper)) {
    hyper <- list()
  }
  if (!is.list(hyper)) {
    stop("'hyper' must be a list.", call. = FALSE)
  }
  if (!is.null(hyper$prec)) {
    stop(
      "Use 'hyper$kappa', not 'hyper$prec', for the LAvM concentration.",
      call. = FALSE
    )
  }

  kappa.setting <- hyper$kappa
  if (is.null(kappa.setting)) {
    kappa.setting <- list()
  }
  if (!is.list(kappa.setting)) {
    stop("'hyper$kappa' must be a list.", call. = FALSE)
  }
  if (!is.null(kappa.setting$log.initial)) {
    stop(
      paste0(
        "'log.initial' is no longer used for LAvM. Use 'initial'; ",
        "it is already interpreted on the log(kappa) scale."
      ),
      call. = FALSE
    )
  }

  prior.name <- if (is.null(kappa.setting$prior)) {
    "pc.vminf"
  } else {
    .canonical_lavm_pc_prior(kappa.setting$prior)
  }

  prior.param <- if (is.null(kappa.setting$param)) {
    c(0.5, 0.5)
  } else {
    as.numeric(kappa.setting$param)
  }
  if (length(prior.param) != 2L || anyNA(prior.param) ||
      any(!is.finite(prior.param))) {
    stop(
      "'hyper$kappa$param' must be the two finite values c(u, alpha).",
      call. = FALSE
    )
  }
  prior.u <- prior.param[1L]
  prior.alpha <- prior.param[2L]
  if (prior.u <= 0 || prior.u >= 1) {
    stop("The PC-prior parameter 'u' must lie strictly between 0 and 1.",
         call. = FALSE)
  }
  if (prior.alpha <= 0 || prior.alpha >= 1) {
    stop(
      "The PC-prior parameter 'alpha' must lie strictly between 0 and 1.",
      call. = FALSE
    )
  }

  initial.theta <- if (is.null(kappa.setting$initial)) {
    6.0
  } else {
    as.numeric(kappa.setting$initial)
  }
  if (length(initial.theta) != 1L || is.na(initial.theta) ||
      !is.finite(initial.theta)) {
    stop(
      "'hyper$kappa$initial' must be one finite value on the log(kappa) scale.",
      call. = FALSE
    )
  }

  fixed.theta <- if (is.null(kappa.setting$fixed)) {
    FALSE
  } else {
    kappa.setting$fixed
  }
  if (!is.logical(fixed.theta) || length(fixed.theta) != 1L ||
      is.na(fixed.theta)) {
    stop("'hyper$kappa$fixed' must be TRUE or FALSE.", call. = FALSE)
  }

  list(
    setting = family.setting,
    link = link,
    link.code = link.code,
    prior = prior.name,
    prior.code = if (prior.name == "pc.vm0") 1L else 0L,
    u = prior.u,
    alpha = prior.alpha,
    initial = initial.theta,
    fixed = isTRUE(fixed.theta)
  )
}

.INLAcircular_shlib_path <- function() {
  loaded <- getLoadedDLLs()
  if ("INLAcircular" %in% names(loaded)) {
    path <- tryCatch(
      loaded[["INLAcircular"]][["path"]],
      error = function(e) ""
    )
    if (is.character(path) && length(path) == 1L && nzchar(path)) {
      return(path)
    }
  }

  path <- system.file(
    "libs",
    paste0("INLAcircular", .Platform$dynlib.ext),
    package = "INLAcircular"
  )
  if (!nzchar(path)) {
    path <- paste0("src/INLAcircular", .Platform$dynlib.ext)
  }
  path
}

.define_lavm_cloglike <- function(info) {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "Package 'INLA' is required. Install the INLA testing version first.",
      call. = FALSE
    )
  }

  INLA::inla.cloglike.define(
    model = "INLAcirc_cloglike_lavm",
    shlib = .INLAcircular_shlib_path(),
    lavm.link = as.numeric(info$link.code),
    lavm.prior = as.numeric(info$prior.code),
    lavm.u = as.numeric(info$u),
    lavm.alpha = as.numeric(info$alpha),
    lavm.initial.theta = as.numeric(info$initial),
    lavm.fixed.theta = as.numeric(info$fixed)
  )
}

.lavm_prior_label <- function(info) {
  sprintf("%s(%g, %g)", info$prior, info$u, info$alpha)
}

.make_lavm_inla_family <- function(family.setting = NULL) {
  info <- .parse_lavm_family_setting(family.setting)
  setting <- info$setting

  # link and hyper are consumed by the LAvM C module. Other ordinary
  # control.family entries are forwarded to INLA's cloglike family.
  extras <- setting[setdiff(
    names(setting),
    c("link", "hyper", "cloglike", "lambda", "u", "alpha")
  )]
  control <- c(
    list(cloglike = .define_lavm_cloglike(info)),
    extras
  )

  list(
    family = "cloglike",
    control = control,
    info = info,
    prior.label = .lavm_prior_label(info)
  )
}

#' Define the custom LAVM likelihood for INLA
#'
#' @param family.setting A list of controls for the LAvM likelihood. Within
#'   `hyper$kappa`, `initial` is the initial value of `log(kappa)`, and
#'   `fixed = TRUE` fixes `log(kappa)` at that value. The two values in
#'   `param = c(u, alpha)` must both lie strictly between 0 and 1. For a single
#'   likelihood, both `list(link = ..., hyper = ...)` and the nested INLA form
#'   `list(list(link = ..., hyper = ...))` are accepted.
#' @details The native likelihood entry point is
#'   `INLAcirc_cloglike_lavm`. Its implementation is kept in the R-independent
#'   `src-cloglike` source tree. Supported concentration priors are
#'   `pc.vm0` and `pc.vminf`.
#' @return An INLA cloglike object.
#' @export
lavm.cloglike <- function(family.setting = NULL) {
  .define_lavm_cloglike(.parse_lavm_family_setting(family.setting))
}

#' Post-process INLA output to rename and exponentiate LAVM parameters
#' @keywords internal
lavm.rename.inla.output <- function(result, meta) {
  if (is.null(result) || is.null(meta) || !("lavm" %in% meta$base_families)) return(result)

  h_names <- rownames(result$summary.hyperpar)
  idx <- grep("Theta[0-9]+ for (INLA\\.Data|cloglike)", h_names)

  if (length(idx) == 0) return(result)

  # Use seq_along to create a clean index (1, 2, ...)
  for (i in seq_along(idx)) {
    k <- idx[i]
    old_name <- h_names[k]

    # If there are multiple LAVM likelihoods, append [1], [2], etc.
    suffix <- if (length(idx) > 1) paste0("[", i, "]") else ""

    new_ext_name <- paste0("kappa for lavm observations", suffix)
    new_int_name <- paste0("Log kappa for lavm observations", suffix)

    if (old_name %in% names(result$marginals.hyperpar)) {
      marg <- result$marginals.hyperpar[[old_name]]
      marg_kappa <- INLA::inla.tmarginal(function(x) exp(x), marg)

      z <- INLA::inla.zmarginal(marg_kappa, silent = TRUE)
      m <- INLA::inla.mmarginal(marg_kappa)

      result$summary.hyperpar[k, "mean"] <- z$mean
      result$summary.hyperpar[k, "sd"] <- z$sd
      result$summary.hyperpar[k, "0.025quant"] <- z$quant0.025
      result$summary.hyperpar[k, "0.5quant"] <- z$quant0.5
      result$summary.hyperpar[k, "0.975quant"] <- z$quant0.975
      result$summary.hyperpar[k, "mode"] <- m

      result$marginals.hyperpar[[old_name]] <- marg_kappa
    }

    rownames(result$summary.hyperpar)[k] <- new_ext_name
    if (old_name %in% names(result$marginals.hyperpar)) {
      names(result$marginals.hyperpar)[names(result$marginals.hyperpar) == old_name] <- new_ext_name
    }

    int_names <- rownames(result$internal.summary.hyperpar)
    int_idx <- which(int_names == old_name)
    if (length(int_idx) > 0) {
      rownames(result$internal.summary.hyperpar)[int_idx] <- new_int_name
    }
    if (old_name %in% names(result$internal.marginals.hyperpar)) {
      names(result$internal.marginals.hyperpar)[names(result$internal.marginals.hyperpar) == old_name] <- new_int_name
    }

    # Bridge the prior name so cc.summary.R can accurately assign it
    if (!is.null(result$inlacc_meta$hyper_priors)) {
      if (old_name %in% names(result$inlacc_meta$hyper_priors)) {
        val <- result$inlacc_meta$hyper_priors[[old_name]]
        result$inlacc_meta$hyper_priors[[new_ext_name]] <- val
      }
    }
  }
  return(result)
}
