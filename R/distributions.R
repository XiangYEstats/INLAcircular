#' @useDynLib INLAcircular, .registration = TRUE
#' @importFrom stats qlogis plogis runif
NULL

# Hidden environment to securely cache grids
.inlacc_cache <- new.env(parent = emptyenv())

#' Generate Spline Grid for Standard Von Mises (mu = 0)
#'
#' @param kappa Numeric scalar. Concentration.
#' @param len Integer. Grid size.
#' @keywords internal
vm.spline.grid <- function(kappa, len = 2048L) {
  s <- 1 / sqrt(kappa)
  if (s < 1/sqrt(15)) {
    x.grid <- seq(-5*s, 5*s, length.out = len)
  } else {
    x.grid <- seq(-pi, pi, length.out = len)
  }
  fx <- dvm(x.grid, mu = 0, kappa = kappa, log = FALSE)
  Fx <- cumsum(fx)
  Fx.norm <- (Fx - Fx[1]) / (Fx[length(Fx)] - Fx[1])
  keep <- (Fx.norm > 0) & (Fx.norm < 1)
  if (sum(keep) < 4) {
    return(data.frame(x = x.grid, logit.Fx = qlogis(punif(seq(0,1,length.out=length(x.grid))))))
  }
  return(data.frame(x = x.grid[keep], logit.Fx = qlogis(Fx.norm[keep])))
}

#' Helper to get or compute grid and memoize it
#' @keywords internal
get_vm_grid <- function(kappa, len = 2048L) {
  cache_name <- paste0("k_", round(kappa, 5), "_l_", len)

  if (!exists(cache_name, envir =.inlacc_cache)) {
    grid <- vm.spline.grid(kappa, len)
    assign(cache_name, grid, envir =.inlacc_cache)
  }
  return(get(cache_name, envir =.inlacc_cache))
}

#' Create Link-Adjusted von Mises (LAvM) Link Object
#'
#' @param type Character string defining the link function. Default is "tan".
#' @return A list of class \code{lavm_link} containing the transform and inverse functions.
#' @export
lavm.link <- function(type = "tan") {
  if (type == "tan") {
    structure(list(
      name = "tan",
      transform = function(x, eta) {
        y <- tan(x / 2)
        u <- y - eta
        z <- 2 * atan(u)
        log.J <- log1p(y^2) - log1p(u^2)
        return(list(z = z, log.J = log.J))
      },
      inverse = function(z, eta) {
        return(2 * atan(tan(z / 2) + eta))
      }
    ), class = "lavm_link")
  } else {
    stop(sprintf("Link '%s' not implemented", type))
  }
}

# ==============================================================================
# 1. Standard von Mises Distribution
# ==============================================================================

#' The von Mises Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the von Mises distribution with mean direction \code{mu} and
#' concentration parameter \code{kappa}.
#'
#' @param x,q vector of quantiles (angles in radians).
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu mean direction of the distribution (in radians).
#' @param kappa concentration parameter (must be strictly positive).
#' @param log logical; if TRUE, probabilities/densities are given as log.
#' @param strategy character; either \code{"circular"} or \code{"linear"} for the CDF integration strategy. If \code{"circular"}, the CDF integration starts from \eqn{\mu - \pi}. If \code{"linear"}, it starts from \eqn{-\pi}.
#' @param len integer; resolution of the spline grid used for approximation.
#'
#' @details
#' The von Mises distribution has the following density function:
#' \deqn{f\left(x \mid \mu, \kappa\right) = \frac{1}{2\pi I_0(\kappa)} \exp\left\{\kappa \cos\left(x - \mu\right)\right\}}
#' where \eqn{I_0\left(\kappa\right)} is the modified Bessel function of the first kind of order 0.
#'
#' @return \code{dvm} gives the density, \code{pvm} gives the distribution function,
#' \code{qvm} gives the quantile function, and \code{rvm} generates random deviates.
#'
#' @name von_mises
#' @rdname von_mises
#'
#' @examples
#' # Generate 10 random deviates from a von Mises distribution
#' # with mean direction \eqd{\pi/2} and concentration 5
#' rvm(10, mu = pi/2, kappa = 5)
#'
#' # Calculate the density at specific points
#' x_grid <- seq(-pi, pi, length.out = 100)
#' densities <- dvm(x_grid, mu = 0, kappa = 2)
#'
#' # Plot the density curve
#' plot(x_grid, densities, type = "l", main = "von Mises Density (kappa = 2)",
#'      xlab = "Angle", ylab = "Density")
#'
#' # Calculate the cumulative probability up to pi/2
#' pvm(pi/2, mu = 0, kappa = 2)
#'
#' # Find the 95th percentile
#' qvm(0.95, mu = 0, kappa = 2)
#'
#' @export
dvm <- function(x, mu, kappa, log = FALSE) {
  args <- data.frame(x = x, mu = mu, kappa = kappa)
  .Call("INLAcirc_C_dvm",
        as.numeric(args$x),
        as.numeric(args$mu),
        as.numeric(args$kappa),
        as.logical(log),
        PACKAGE = "INLAcircular")
}

#' @rdname von_mises
#' @export
pvm <- function(q, mu, kappa, strategy = "circular", log = FALSE, len = 2048L) {
  args <- data.frame(q = q, mu = mu, kappa = kappa)
  res <- numeric(nrow(args))
  strat_int <- ifelse(strategy == "circular", 0L, 1L)

  for (k in unique(args$kappa)) {
    idx <- which(args$kappa == k)
    grid <- get_vm_grid(k, len)

    res[idx] <-.Call("INLAcirc_C_pvm",
                     as.numeric(args$q[idx]),
                     as.numeric(args$mu[idx]),
                     as.integer(strat_int),
                     as.logical(log),
                     as.numeric(grid$x),
                     as.numeric(grid$logit.Fx),
                     PACKAGE = "INLAcircular")
  }
  return(res)
}

#' @rdname von_mises
#' @export
qvm <- function(p, mu, kappa, len = 2048L) {
  args <- data.frame(p = p, mu = mu, kappa = kappa)
  res <- numeric(nrow(args))

  for (k in unique(args$kappa)) {
    idx <- which(args$kappa == k)
    grid <- get_vm_grid(k, len)

    res[idx] <-.Call("INLAcirc_C_qvm", as.numeric(args$p[idx]), as.numeric(args$mu[idx]),
                     as.numeric(grid$x), as.numeric(grid$logit.Fx),
                     PACKAGE = "INLAcircular")
  }
  return(res)
}

#' @rdname von_mises
#' @export
rvm <- function(n, mu, kappa, len = 2048L) {
  if (length(mu) == 1) mu <- rep(mu, n)
  if (length(kappa) == 1) kappa <- rep(kappa, n)

  res <- numeric(n)
  for (k in unique(kappa)) {
    idx <- which(kappa == k)
    n_idx <- length(idx)
    grid <- get_vm_grid(k, len)

    res[idx] <-.Call("INLAcirc_C_rvm", as.integer(n_idx), as.numeric(mu[idx]),
                     as.numeric(grid$x), as.numeric(grid$logit.Fx),
                     PACKAGE = "INLAcircular")
  }
  return(res)
}


# ==============================================================================
# 2. PC Prior for the von Mises Concentration
# ==============================================================================

.pc_vm_check_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}

.pc_vm_recycle <- function(value, lambda, value_name) {
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

.pc_vm_sample_size <- function(n) {
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
#' The density and CDF are
#' \deqn{\pi_0(\kappa)=\lambda\exp\{-\lambda d_0(\kappa)\}
#'   |d_0'(\kappa)|,}
#' \deqn{F_0(\kappa)=1-\exp\{-\lambda d_0(\kappa)\}.}
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
  log <- .pc_vm_check_flag(log, "log")
  args <- .pc_vm_recycle(kappa, lambda, "kappa")
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
  log.p <- .pc_vm_check_flag(log.p, "log.p")
  args <- .pc_vm_recycle(q, lambda, "q")
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
  log.p <- .pc_vm_check_flag(log.p, "log.p")
  args <- .pc_vm_recycle(p, lambda, "p")
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
  n <- .pc_vm_sample_size(n)
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
#' The continuous density and CDF are
#' \deqn{\pi_\infty(\kappa)=
#'   \lambda\exp\{-\lambda d_\infty(\kappa)\}
#'   |d_\infty'(\kappa)|,}
#' \deqn{F_\infty(\kappa)=\exp\{-\lambda d_\infty(\kappa)\}.}
#'
#' Because \eqn{d_\infty(0)=1}, the CDF contains boundary mass
#' \eqn{\exp(-\lambda)} at \eqn{\kappa=0}. \code{dpc.vminf} returns the
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
  log <- .pc_vm_check_flag(log, "log")
  args <- .pc_vm_recycle(kappa, lambda, "kappa")
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
  log.p <- .pc_vm_check_flag(log.p, "log.p")
  args <- .pc_vm_recycle(q, lambda, "q")
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
  log.p <- .pc_vm_check_flag(log.p, "log.p")
  args <- .pc_vm_recycle(p, lambda, "p")
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
  n <- .pc_vm_sample_size(n)
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


# ==============================================================================
# 3. Link-Adjusted von Mises (LAvM) Distribution
# ==============================================================================

#' The Link-Adjusted von Mises (LAvM) Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the Link-Adjusted von Mises (LAvM) distribution with linear predictor \code{eta}
#' and concentration parameter \code{kappa}.
#'
#' @param x,q vector of quantiles (angles in radians).
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param eta linear predictor vector.
#' @param kappa concentration parameter (must be strictly positive).
#' @param log logical; if TRUE, probabilities/densities are given as log.
#' @param link_obj optional link object generated by \code{lavm.link()}.
#' @param len integer; resolution of the spline grid used for approximation.
#'
#' @details
#' The default transformation from the standard circular space to the linear real line
#' is given by the \emph{inverse tangent link},
#' \deqn{g\left(\cdot\right) := 2 \arctan\left(\cdot\right) \in \left[-\pi, \pi \right),}
#' resulting in:
#' \deqn{z = 2 \arctan\left( \tan\left( \frac{x}{2} \right) - \eta \right)}
#'
#' The corresponding density function for the Link-Adjusted von Mises (LAvM) distribution is:
#' \deqn{
#'   \pi_{\operatorname{LAvM}}\left(x \mid \eta, \kappa\right) =
#'   \frac{\exp\left\lbrace \kappa \cos\left( 2 \arctan \left( \tan \left(
#'   \frac{x}{2} \right) - \eta \right) \right) \right\rbrace}{2\pi I_{0}
#'   \left(\kappa\right)\left(1 + \eta^{2} - \eta \sin\left(x\right) -
#'   \eta^{2}\sin^{2}\left(\frac{x}{2}\right)\right)},
#'   \quad x\in\left(-\pi,\pi\right)
#' }
#'
#' @return \code{dlavm} gives the density, \code{plavm} gives the distribution function,
#' \code{qlavm} gives the quantile function, and \code{rlavm} generates random deviates.
#'
#' @name lavm
#' @rdname lavm
#'
#' @examples
#' # Generate 1000 random deviates from the LAvM distribution
#' # with linear predictor eta = 1.5 and concentration kappa = 3
#' rlavm(10, eta = 1.5, kappa = 3)
#'
#' # Calculate the density at specific points
#' x_grid <- seq(-pi, pi, length.out = 100)
#' densities <- dlavm(x_grid, eta = 1.5, kappa = 3)
#'
#' # Plot the density curve
#' plot(x_grid, densities, type = "l", main = "LAvM Density (eta = 1.5, kappa = 3)",
#'      xlab = "Angle", ylab = "Density")
#'
#' # Calculate the cumulative probability up to pi/2
#' plavm(pi/2, eta = 1.5, kappa = 3)
#'
#' # Find the median (50th percentile)
#' qlavm(0.5, eta = 1.5, kappa = 3)
#'
#' @export
dlavm <- function(x, eta, kappa, log = FALSE) {
  args <- data.frame(x = x, eta = eta, kappa = kappa)
  .Call("INLAcirc_C_dlavm",
        as.numeric(args$x),
        as.numeric(args$eta),
        as.numeric(args$kappa),
        as.logical(log),
        PACKAGE = "INLAcircular")
}

#' @rdname lavm
#' @export
plavm <- function(q, eta, kappa, link_obj = NULL, log = FALSE, len = 2048L) {
  args <- data.frame(q = q, eta = eta, kappa = kappa)
  res <- numeric(nrow(args))

  for (k in unique(args$kappa)) {
    idx <- which(args$kappa == k)
    grid <- get_vm_grid(k, len)

    res[idx] <-.Call("INLAcirc_C_plavm", as.numeric(args$q[idx]), as.numeric(args$eta[idx]),
                     as.logical(log), as.numeric(grid$x), as.numeric(grid$logit.Fx),
                     PACKAGE = "INLAcircular")
  }
  return(res)
}

#' @rdname lavm
#' @export
qlavm <- function(p, eta, kappa, link_obj = NULL, len = 2048L) {
  args <- data.frame(p = p, eta = eta, kappa = kappa)
  res <- numeric(nrow(args))

  for (k in unique(args$kappa)) {
    idx <- which(args$kappa == k)
    grid <- get_vm_grid(k, len)

    res[idx] <-.Call("INLAcirc_C_qlavm", as.numeric(args$p[idx]), as.numeric(args$eta[idx]),
                     as.numeric(grid$x), as.numeric(grid$logit.Fx),
                     PACKAGE = "INLAcircular")
  }
  return(res)
}

#' @rdname lavm
#' @export
rlavm <- function(n, eta, kappa, link_obj = NULL, len = 2048L) {
  if (length(eta) == 1) eta <- rep(eta, n)
  if (length(kappa) == 1) kappa <- rep(kappa, n)

  res <- numeric(n)
  for (k in unique(kappa)) {
    idx <- which(kappa == k)
    n_idx <- length(idx)
    grid <- get_vm_grid(k, len)

    res[idx] <-.Call("INLAcirc_C_rlavm", as.integer(n_idx), as.numeric(eta[idx]),
                     as.numeric(grid$x), as.numeric(grid$logit.Fx),
                     PACKAGE = "INLAcircular")
  }
  return(res)
}
