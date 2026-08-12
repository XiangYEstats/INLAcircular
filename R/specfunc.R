#' Modified Bessel Function of the First Kind
#'
#' Evaluates the modified Bessel function of the first kind using high-precision
#' Chebyshev approximations.
#'
#' @param x numeric vector of values.
#' @param nu numeric scalar indicating the order of the Bessel function.
#' @param expon.scaled logical; if TRUE, the results are exponentially scaled
#'        by \code{exp(-abs(x))}. Default is FALSE.
#' @return A numeric vector of calculated Bessel values.
#' @export
bessel_i <- function(x, nu, expon.scaled = FALSE) {
  if (length(nu) != 1) stop("'nu' must be a scalar.")

  .Call("INLAcirc_C_bessel_i",
        as.numeric(x),
        as.numeric(nu),
        as.logical(expon.scaled),
        PACKAGE = "INLAcircular")
}
