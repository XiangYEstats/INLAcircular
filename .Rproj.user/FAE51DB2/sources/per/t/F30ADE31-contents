#' Define the custom LAVM likelihood for INLA
#'
#' @param family.setting A list of controls for the LAVM likelihood.
#' @return An INLA cloglike object.
#' @export
lavm.cloglike <- function(family.setting = NULL) {

  link <- if (!is.null(family.setting$link)) as.character(family.setting$link) else "inverse.tangent"

  link.code <- switch(link,
                      "inverse.tangent" = 0L,
                      "scaled.logit"    = 1L,
                      "scaled.probit"   = 2L)

  if (is.null(link.code)) stop("Invalid link function.")

  prior.u <- 0.5
  prior.alpha <- 0.5
  initial.theta <- 6.0
  prior.name <- "pc.vm.inf"

  if (!is.null(family.setting$hyper) && !is.null(family.setting$hyper$kappa)) {
    kappa_set <- family.setting$hyper$kappa

    if (!is.null(kappa_set$prior)) {
      prior.name <- as.character(kappa_set$prior)
    }

    if (!is.null(kappa_set$log.initial)) {
      initial.theta <- as.numeric(kappa_set$log.initial)
    } else if (!is.null(kappa_set$initial)) {
      initial.theta <- log(as.numeric(kappa_set$initial))
    }

    if (!is.null(kappa_set$param) && length(kappa_set$param) >= 2) {
      prior.u <- as.numeric(kappa_set$param[1])
      prior.alpha <- as.numeric(kappa_set$param[2])
    }
  }

  prior.code <- if (prior.name == "pc.vm.0") 1L else 0L

  shlib.path <- system.file("libs", paste0("INLAcircular", .Platform$dynlib.ext), package = "INLAcircular")
  if (shlib.path == "") {
    shlib.path <- paste0("src/INLAcircular", .Platform$dynlib.ext)
  }

  INLA::inla.cloglike.define(
    model = "inla_cloglike_lavm",
    shlib = shlib.path,
    lavm.link = as.numeric(link.code),
    lavm.prior = as.numeric(prior.code),
    lavm.u = as.numeric(prior.u),
    lavm.alpha = as.numeric(prior.alpha),
    lavm.initial.theta = as.numeric(initial.theta)
  )
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
