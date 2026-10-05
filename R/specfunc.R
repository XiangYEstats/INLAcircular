#' Modified Bessel Function of the First Kind
#'
#' Evaluates the modified Bessel function of the first kind using native C
#' routines adapted from the GNU Scientific Library (GSL). R supplies only the
#' interface; argument checking and numerical computation take place in C.
#'
#' @param x Numeric vector. Negative values are supported for integer orders.
#'   Missing values are propagated.
#' @param nu One finite, nonnegative numeric value indicating the order.
#'   Fractional orders are supported for nonnegative `x`.
#' @param expon.scaled Logical scalar; if `TRUE`, return
#'   \eqn{I_\nu(x)\exp(-|x|)}. Defaults to `FALSE`.
#' @details
#' Negative `x` with a fractional order produces `NaN` and a warning. For an
#' integer order \eqn{n}, \eqn{I_n(-x) = (-1)^n I_n(x)}.
#'
#' At zero, the result is one for order zero and zero for positive orders.
#' At positive infinity, the result is zero when scaled and infinity otherwise.
#' Unscaled results can overflow for large arguments; use `expon.scaled = TRUE`
#' when the scaled value is needed.
#' For very large orders, cancellation can reduce the accuracy of an unscaled
#' result. A warning is issued when the estimated loss is significant.
#'
#' @return A numeric vector of calculated Bessel values.
#' @examples
#' bessel_i(1, nu = 0.5)
#' bessel_i(c(0, 1, 5), nu = 0)
#' bessel_i(1000, nu = 0, expon.scaled = TRUE)
#' @export
bessel_i <- function(x, nu, expon.scaled = FALSE) {
  .Call("C_bessel_i", x, nu, expon.scaled, PACKAGE = "INLAcircular")
}
