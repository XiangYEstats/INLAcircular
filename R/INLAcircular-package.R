#' INLAcircular: Bayesian Joint Circular Regression with INLA
#'
#' \pkg{INLAcircular} supplies compiled circular probability distributions,
#' penalized-complexity priors for their concentration parameters, and a
#' formula interface for Bayesian circular and joint circular-linear models
#' fitted with Integrated Nested Laplace Approximation (INLA).
#'
#' The main entry points are:
#' \itemize{
#'   \item \code{\link{dvm}} and the other \code{d/p/q/r} functions for the
#'     von Mises and link-adjusted von Mises distributions;
#'   \item \code{\link{dcardioid}} and \code{\link{dwrappedcauchy}} for
#'     cardioid and wrapped Cauchy densities;
#'   \item \code{\link{dpc.vm0}}, \code{\link{dpc.vminf}},
#'     \code{\link{dpc.card0}}, \code{\link{dpc.card}}, and
#'     \code{\link{dpc.wc}} families of PC-prior functions;
#'   \item \code{\link{likelihood}}, \code{\link{inlacc}}, and
#'     \code{\link{inla}} for model construction and fitting; and
#'   \item \code{\link{LKJcc}} and \code{\link{index}} for multivariate and
#'     structured latent effects.
#' }
#'
#' See the package guide, \code{vignette("INLAcircular-guide")}, for a
#' complete workflow and function index.
#'
#' @docType package
#' @name INLAcircular-package
#' @aliases INLAcircular
#' @keywords internal
#' @importFrom stats as.formula punif terms
"_PACKAGE"
