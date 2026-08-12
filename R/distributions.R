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
# 2. Link-Adjusted von Mises (LAvM) Distribution
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
