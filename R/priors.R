# ==============================================================================
# PC Priors for Circular Concentration Parameters
# ==============================================================================

.pc_check_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}

.pc_recycle <- function(value, lambda, value_name) {
  value <- as.numeric(value)
  lambda <- as.numeric(lambda)

  if (length(lambda) == 0L || anyNA(lambda) ||
      any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("'lambda' must contain positive, finite values.", call. = FALSE)
  }
  if (length(value) == 0L) {
    return(list(value = numeric(), lambda = numeric()))
  }

  output_length <- max(length(value), length(lambda))
  if (output_length %% length(value) != 0L ||
      output_length %% length(lambda) != 0L) {
    warning(sprintf(
      "longer object length is not a multiple of shorter object length in '%s' and 'lambda'",
      value_name
    ), call. = FALSE)
  }

  list(
    value = rep_len(value, output_length),
    lambda = rep_len(lambda, output_length)
  )
}

.pc_sample_size <- function(n) {
  if (length(n) > 1L) {
    return(length(n))
  }
  if (length(n) != 1L || is.na(n) || !is.finite(n) || n < 0) {
    stop("'n' must be a non-negative finite number.", call. = FALSE)
  }

  n <- as.integer(n)
  if (is.na(n)) {
    stop("'n' is too large.", call. = FALSE)
  }
  n
}

#' PC Prior for von Mises Concentration: Uniform Base
#'
#' Density, distribution function, quantile function, and random generation
#' for the penalized-complexity prior on the von Mises concentration
#' parameter with the circular uniform distribution (\eqn{\kappa=0}) as its
#' base model.
#'
#' @param kappa,q Vector of non-negative concentration values.
#' @param p Vector of probabilities.
#' @param n Number of observations. If \code{length(n) > 1}, its length is
#'   used.
#' @param lambda Positive rate parameter.
#' @param log Logical; if \code{TRUE}, return the log-density.
#' @param log.p Logical; if \code{TRUE}, probabilities are supplied or
#'   returned on the log scale.
#'
#' @details
#' The PC distance from the uniform base model is
#' \deqn{d_0(\kappa) =
#'   \sqrt{\kappa I_1(\kappa)/I_0(\kappa)-\log I_0(\kappa)}.}
#' Write \eqn{A(\kappa)=I_1(\kappa)/I_0(\kappa)} and
#' \eqn{A'(\kappa)=1-A(\kappa)^2-A(\kappa)/\kappa}, with the continuous
#' limit \eqn{A'(0)=1/2}. The density and CDF are
#' \deqn{\pi_0(\kappa)=\lambda\exp\{-\lambda d_0(\kappa)\}
#'   \frac{\kappa A'(\kappa)}{2d_0(\kappa)},}
#' \deqn{F_0(\kappa)=1-\exp\{-\lambda d_0(\kappa)\}.}
#' The density has the finite endpoint limit \eqn{\pi_0(0)=\lambda/2};
#' there is no atom at zero.
#'
#' All Bessel evaluations, small-\eqn{\kappa} density expansions,
#' large-\eqn{\kappa} asymptotic expansions, CDF calculations, and numerical
#' inversion are performed by the package's compiled C code.
#'
#' @return \code{dpc.vm0} gives the density, \code{ppc.vm0} gives the
#'   distribution function, \code{qpc.vm0} gives the quantile function, and
#'   \code{rpc.vm0} generates random deviates.
#'
#' @name pc_vm0
#' @rdname pc_vm0
#' @aliases pc.vm0
#'
#' @examples
#' dpc.vm0(kappa = c(0, 1, 5), lambda = 2)
#' ppc.vm0(q = 5, lambda = 2)
#' qpc.vm0(p = 0.5, lambda = 2)
#' rpc.vm0(n = 10, lambda = 2)
#'
#' @export
dpc.vm0 <- function(kappa, lambda, log = FALSE) {
  log <- .pc_check_flag(log, "log")
  args <- .pc_recycle(kappa, lambda, "kappa")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_dpc_vm0",
    args$value,
    args$lambda,
    log,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_vm0
#' @export
ppc.vm0 <- function(q, lambda, log.p = FALSE) {
  log.p <- .pc_check_flag(log.p, "log.p")
  args <- .pc_recycle(q, lambda, "q")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_ppc_vm0",
    args$value,
    args$lambda,
    log.p,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_vm0
#' @export
qpc.vm0 <- function(p, lambda, log.p = FALSE) {
  log.p <- .pc_check_flag(log.p, "log.p")
  args <- .pc_recycle(p, lambda, "p")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_qpc_vm0",
    args$value,
    args$lambda,
    log.p,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_vm0
#' @export
rpc.vm0 <- function(n, lambda) {
  n <- .pc_sample_size(n)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 0L || anyNA(lambda) ||
      any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("'lambda' must contain positive, finite values.", call. = FALSE)
  }
  if (n == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_rpc_vm0",
    n,
    rep_len(lambda, n),
    PACKAGE = "INLAcircular"
  )
}

#' PC Prior for von Mises Concentration: Point-Mass Base
#'
#' Density, distribution function, quantile function, and random generation
#' for the penalized-complexity prior on the von Mises concentration
#' parameter with the point-mass limit (\eqn{\kappa\to\infty}) as its base
#' model.
#'
#' @inheritParams pc_vm0
#'
#' @details
#' The PC distance from the point-mass base model is
#' \deqn{d_\infty(\kappa)=
#'   \sqrt{1-I_1(\kappa)/I_0(\kappa)}.}
#' Writing \eqn{A(\kappa)=I_1(\kappa)/I_0(\kappa)} and
#' \eqn{A'(\kappa)=1-A(\kappa)^2-A(\kappa)/\kappa}, the continuous density
#' and full CDF are
#' \deqn{\pi_\infty^{\mathrm{cont}}(\kappa)=
#'   \lambda\exp\{-\lambda d_\infty(\kappa)\}
#'   \frac{A'(\kappa)}{2d_\infty(\kappa)},}
#' \deqn{F_\infty(\kappa)=\exp\{-\lambda d_\infty(\kappa)\}.}
#'
#' Because \eqn{d_\infty(0)=1}, the CDF contains boundary mass
#' \eqn{\exp(-\lambda)} at \eqn{\kappa=0}. The continuous density has
#' right-hand limit \eqn{\lambda\exp(-\lambda)/4} there.
#' \code{dpc.vminf} returns the
#' continuous density component; \code{ppc.vminf}, \code{qpc.vminf}, and
#' \code{rpc.vminf} include the boundary mass. All numerical calculations,
#' including the large-\eqn{\kappa} asymptotic branch and quantile inversion,
#' are performed by compiled C code.
#'
#' @return \code{dpc.vminf} gives the continuous density,
#'   \code{ppc.vminf} gives the distribution function,
#'   \code{qpc.vminf} gives the quantile function, and \code{rpc.vminf}
#'   generates random deviates.
#'
#' @name pc_vminf
#' @rdname pc_vminf
#' @aliases pc.vminf
#'
#' @examples
#' dpc.vminf(kappa = c(0, 1, 5), lambda = 2)
#' ppc.vminf(q = 5, lambda = 2)
#' qpc.vminf(p = 0.5, lambda = 2)
#' rpc.vminf(n = 10, lambda = 2)
#'
#' @export
dpc.vminf <- function(kappa, lambda, log = FALSE) {
  log <- .pc_check_flag(log, "log")
  args <- .pc_recycle(kappa, lambda, "kappa")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_dpc_vminf",
    args$value,
    args$lambda,
    log,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_vminf
#' @export
ppc.vminf <- function(q, lambda, log.p = FALSE) {
  log.p <- .pc_check_flag(log.p, "log.p")
  args <- .pc_recycle(q, lambda, "q")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_ppc_vminf",
    args$value,
    args$lambda,
    log.p,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_vminf
#' @export
qpc.vminf <- function(p, lambda, log.p = FALSE) {
  log.p <- .pc_check_flag(log.p, "log.p")
  args <- .pc_recycle(p, lambda, "p")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_qpc_vminf",
    args$value,
    args$lambda,
    log.p,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_vminf
#' @export
rpc.vminf <- function(n, lambda) {
  n <- .pc_sample_size(n)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 0L || anyNA(lambda) ||
      any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("'lambda' must contain positive, finite values.", call. = FALSE)
  }
  if (n == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_rpc_vminf",
    n,
    rep_len(lambda, n),
    PACKAGE = "INLAcircular"
  )
}


#' PC Prior for Cardioid Concentration: Uniform Base
#'
#' Density, distribution function, quantile function, and random generation
#' for the penalized-complexity prior on cardioid concentration with the
#' circular uniform distribution (\eqn{\kappa=0}) as its base model.
#'
#' @param kappa,q Vectors of cardioid concentration values.
#' @param p Vector of probabilities.
#' @param n Number of observations. If \code{length(n) > 1}, its length is
#'   used.
#' @param lambda Positive rate parameter on the PC-distance scale.
#' @param log Logical scalar; if \code{TRUE}, return log-densities.
#' @param log.p Logical scalar; if \code{TRUE}, probabilities are supplied or
#'   returned on the log scale.
#'
#' @details
#' For \eqn{0\leq\kappa\leq 1/2}, let
#' \deqn{s(\kappa)=\sqrt{1-4\kappa^2}}
#' and
#' \deqn{d_0(\kappa)=\sqrt{1-s(\kappa)+
#'   \log\{(1+s(\kappa))/2\}}.}
#' The largest attainable distance is
#' \eqn{D_0=\sqrt{1-\log 2}}. Consequently, the exponential PC prior on
#' distance is truncated and normalized on \eqn{[0,D_0]}:
#' \deqn{\pi_0(\kappa)=
#'   \frac{\lambda e^{-\lambda d_0(\kappa)}}
#'   {1-e^{-\lambda D_0}}
#'   \frac{2\kappa}{(1+s(\kappa))d_0(\kappa)},}
#' \deqn{F_0(\kappa)=
#'   \frac{1-e^{-\lambda d_0(\kappa)}}
#'   {1-e^{-\lambda D_0}}.}
#'
#' The endpoint densities are
#' \eqn{\lambda/\{1-\exp(-\lambda D_0)\}} at zero and
#' \eqn{\lambda\exp(-\lambda D_0)/
#' \{D_0[1-\exp(-\lambda D_0)]\}} at \eqn{\kappa=1/2}.
#' Algebraic reductions, small-\code{kappa} series, \code{log1p}/\code{expm1}
#' calculations, numerical inversion, and random generation are all performed
#' in compiled C code.
#'
#' @return \code{dpc.card0} gives the density, \code{ppc.card0} gives the
#'   distribution function, \code{qpc.card0} gives the quantile function, and
#'   \code{rpc.card0} generates random deviates.
#'
#' @name pc_card0
#' @aliases pc.card0
#'
#' @examples
#' dpc.card0(kappa = c(0, 0.25, 0.5), lambda = 2)
#' ppc.card0(q = 0.25, lambda = 2)
#' qpc.card0(p = 0.5, lambda = 2)
#' rpc.card0(n = 10, lambda = 2)
#'
#' @export
dpc.card0 <- function(kappa, lambda, log = FALSE) {
  log <- .pc_check_flag(log, "log")
  args <- .pc_recycle(kappa, lambda, "kappa")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_dpc_card0",
    args$value,
    args$lambda,
    log,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_card0
#' @export
ppc.card0 <- function(q, lambda, log.p = FALSE) {
  log.p <- .pc_check_flag(log.p, "log.p")
  args <- .pc_recycle(q, lambda, "q")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_ppc_card0",
    args$value,
    args$lambda,
    log.p,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_card0
#' @export
qpc.card0 <- function(p, lambda, log.p = FALSE) {
  log.p <- .pc_check_flag(log.p, "log.p")
  args <- .pc_recycle(p, lambda, "p")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_qpc_card0",
    args$value,
    args$lambda,
    log.p,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_card0
#' @export
rpc.card0 <- function(n, lambda) {
  n <- .pc_sample_size(n)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 0L || anyNA(lambda) ||
      any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("'lambda' must contain positive, finite values.", call. = FALSE)
  }
  if (n == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_rpc_card0",
    n,
    rep_len(lambda, n),
    PACKAGE = "INLAcircular"
  )
}

#' PC Prior for Cardioid Concentration: Cardioid Base
#'
#' Density, distribution function, quantile function, and random generation
#' for the penalized-complexity prior on cardioid concentration with the
#' maximally concentrated cardioid model (\eqn{\kappa=1/2}) as its base.
#'
#' @inheritParams pc_card0
#'
#' @details
#' For \eqn{0\leq\kappa\leq 1/2}, let
#' \deqn{s(\kappa)=\sqrt{1-4\kappa^2}}
#' and
#' \deqn{d_{1/2}(\kappa)=\sqrt{1-2\kappa-s(\kappa)+
#'   \log(1+s(\kappa))}.}
#' This distance decreases from \eqn{D_{1/2}=\sqrt{\log 2}} at zero to zero
#' at the base model. The continuous density and CDF are
#' \deqn{\pi_{1/2}(\kappa)=\lambda e^{-\lambda d_{1/2}(\kappa)}
#'   \frac{1+s(\kappa)-2\kappa}
#'   {(1+s(\kappa))d_{1/2}(\kappa)},}
#' \deqn{F_{1/2}(\kappa)=e^{-\lambda d_{1/2}(\kappa)}.}
#'
#' Since \eqn{F_{1/2}(0)=e^{-\lambda\sqrt{\log 2}}}, the distribution has a
#' boundary atom of that size at \code{kappa = 0}. \code{dpc.card} returns the
#' continuous density only; \code{ppc.card}, \code{qpc.card}, and
#' \code{rpc.card} include the atom. The continuous density has an integrable
#' singularity at \code{kappa = 0.5}, where it is returned as \code{Inf}.
#' At zero, the continuous density is
#' \eqn{\lambda\exp\{-\lambda\sqrt{\log 2}\}/\sqrt{\log 2}}; this is
#' separate from the atom probability.
#'
#' All endpoint series, log-scale calculations, inversion, and random
#' generation are implemented in compiled C code.
#'
#' @return \code{dpc.card} gives the continuous density,
#'   \code{ppc.card} gives the distribution function, \code{qpc.card} gives
#'   the quantile function, and \code{rpc.card} generates random deviates.
#'
#' @name pc_card
#' @aliases pc.card
#'
#' @examples
#' dpc.card(kappa = c(0, 0.25, 0.49), lambda = 2)
#' ppc.card(q = 0.25, lambda = 2)
#' qpc.card(p = 0.5, lambda = 2)
#' rpc.card(n = 10, lambda = 2)
#'
#' @export
dpc.card <- function(kappa, lambda, log = FALSE) {
  log <- .pc_check_flag(log, "log")
  args <- .pc_recycle(kappa, lambda, "kappa")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_dpc_card",
    args$value,
    args$lambda,
    log,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_card
#' @export
ppc.card <- function(q, lambda, log.p = FALSE) {
  log.p <- .pc_check_flag(log.p, "log.p")
  args <- .pc_recycle(q, lambda, "q")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_ppc_card",
    args$value,
    args$lambda,
    log.p,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_card
#' @export
qpc.card <- function(p, lambda, log.p = FALSE) {
  log.p <- .pc_check_flag(log.p, "log.p")
  args <- .pc_recycle(p, lambda, "p")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_qpc_card",
    args$value,
    args$lambda,
    log.p,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_card
#' @export
rpc.card <- function(n, lambda) {
  n <- .pc_sample_size(n)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 0L || anyNA(lambda) ||
      any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("'lambda' must contain positive, finite values.", call. = FALSE)
  }
  if (n == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_rpc_card",
    n,
    rep_len(lambda, n),
    PACKAGE = "INLAcircular"
  )
}

#' PC Prior for Wrapped Cauchy Concentration
#'
#' Density, distribution function, quantile function, and random generation
#' for the penalized-complexity prior on wrapped Cauchy concentration with the
#' circular uniform distribution (\eqn{\kappa=0}) as its base model.
#'
#' @param kappa,q Vectors of wrapped Cauchy concentration values.
#' @inheritParams pc_card0
#'
#' @details
#' For \eqn{0\leq\kappa<1}, the PC distance is
#' \deqn{d_{\mathrm{wc}}(\kappa)=
#'   \sqrt{-\log(1-\kappa^2)}.}
#' It maps the parameter support to the full non-negative distance axis, so
#' the density and CDF require no truncation:
#' \deqn{\pi_{\mathrm{wc}}(\kappa)=
#'   \lambda e^{-\lambda d_{\mathrm{wc}}(\kappa)}
#'   \frac{\kappa}
#'   {(1-\kappa^2)d_{\mathrm{wc}}(\kappa)},}
#' \deqn{F_{\mathrm{wc}}(\kappa)=
#'   1-e^{-\lambda d_{\mathrm{wc}}(\kappa)}.}
#'
#' The density limit at zero is \eqn{\lambda}; it diverges as
#' \code{kappa} approaches one and is returned as \code{Inf} at that limiting
#' endpoint. The analytic quantile is evaluated with stable
#' \code{log1p}/\code{expm1} calculations. Interior-probability quantiles
#' (\eqn{0 < p < 1}) that would round to exactly one are kept at the largest
#' representable value below one.
#'
#' All density, CDF, quantile, and random-generation calculations are
#' performed in compiled C code.
#'
#' @return \code{dpc.wc} gives the density, \code{ppc.wc} gives the
#'   distribution function, \code{qpc.wc} gives the quantile function, and
#'   \code{rpc.wc} generates random deviates.
#'
#' @name pc_wc
#' @aliases pc.wc
#'
#' @examples
#' dpc.wc(kappa = c(0, 0.5, 0.9), lambda = 2)
#' ppc.wc(q = 0.5, lambda = 2)
#' qpc.wc(p = 0.5, lambda = 2)
#' rpc.wc(n = 10, lambda = 2)
#'
#' @export
dpc.wc <- function(kappa, lambda, log = FALSE) {
  log <- .pc_check_flag(log, "log")
  args <- .pc_recycle(kappa, lambda, "kappa")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_dpc_wc",
    args$value,
    args$lambda,
    log,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_wc
#' @export
ppc.wc <- function(q, lambda, log.p = FALSE) {
  log.p <- .pc_check_flag(log.p, "log.p")
  args <- .pc_recycle(q, lambda, "q")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_ppc_wc",
    args$value,
    args$lambda,
    log.p,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_wc
#' @export
qpc.wc <- function(p, lambda, log.p = FALSE) {
  log.p <- .pc_check_flag(log.p, "log.p")
  args <- .pc_recycle(p, lambda, "p")
  if (length(args$value) == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_qpc_wc",
    args$value,
    args$lambda,
    log.p,
    PACKAGE = "INLAcircular"
  )
}

#' @rdname pc_wc
#' @export
rpc.wc <- function(n, lambda) {
  n <- .pc_sample_size(n)
  lambda <- as.numeric(lambda)
  if (length(lambda) == 0L || anyNA(lambda) ||
      any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("'lambda' must contain positive, finite values.", call. = FALSE)
  }
  if (n == 0L) {
    return(numeric())
  }

  .Call(
    "INLAcirc_C_rpc_wc",
    n,
    rep_len(lambda, n),
    PACKAGE = "INLAcircular"
  )
}
