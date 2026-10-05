suppressPackageStartupMessages(library(INLAcircular))

assert_relative <- function(actual, expected, tolerance = 1e-10) {
  stopifnot(
    length(actual) == length(expected),
    all(is.finite(actual)),
    all(is.finite(expected)),
    all(abs(actual - expected) <=
          tolerance * pmax(abs(expected), .Machine$double.xmin))
  )
}

assert_error <- function(expr) {
  result <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(result, "error"))
}

# I_(1/2)(x) has an elementary closed form, independent of a Bessel library.
half_x <- c(1e-12, 1e-5, 0.1, 1, 10, 100, 700)
assert_relative(
  bessel_i(half_x, 0.5),
  sqrt(2 / (pi * half_x)) * sinh(half_x)
)
assert_relative(
  bessel_i(half_x, 0.5, expon.scaled = TRUE),
  -expm1(-2 * half_x) / sqrt(2 * pi * half_x)
)

# Scaling must remain useful when the unscaled function overflows.
large_x <- c(1000, 1e5, 1e8)
assert_relative(
  bessel_i(large_x, 0.5, expon.scaled = TRUE),
  1 / sqrt(2 * pi * large_x)
)
stopifnot(all(suppressWarnings(bessel_i(large_x, 0.5)) == Inf))

# A scaled intermediate can underflow even when the unscaled answer is
# representable. Reference calculated with 70-digit arithmetic (mpmath).
assert_relative(bessel_i(1000, 1500), 388.89959418735873)
stopifnot(identical(bessel_i(1000, 1500, TRUE), 0))
assert_relative(
  bessel_i(1e308, 0, TRUE),
  (1 / sqrt(2 * pi)) / sqrt(1e308)
)

# Cross-check both integral and fractional orders against R's independent
# implementation, with arguments well inside its supported numerical range.
grid_x <- c(1e-4, 0.01, 0.2, 0.8, 2, 5, 20, 100, 500)
for (order in c(0, 0.25, 0.5, 1, 1.5, 2, 5, 10, 25.5, 50)) {
  for (scaled in c(FALSE, TRUE)) {
    assert_relative(
      bessel_i(grid_x, order, expon.scaled = scaled),
      besselI(grid_x, order, expon.scaled = scaled)
    )
  }
}

# Neighboring orders obey I_(nu-1) - I_(nu+1) = (2 nu / x) I_nu.
# Exponential scaling cancels from the same identity.
recurrence_x <- c(0.1, 1, 3, 10, 30)
for (order in c(1, 1.5, 3, 5.25)) {
  assert_relative(
    bessel_i(recurrence_x, order - 1, TRUE) -
      bessel_i(recurrence_x, order + 1, TRUE),
    (2 * order / recurrence_x) * bessel_i(recurrence_x, order, TRUE),
    tolerance = 1e-9
  )
}

# Integer-order parity applies on the negative real axis, including limits.
parity_x <- c(0.01, 0.5, 2, 10, 100)
for (order in 0:5) {
  for (scaled in c(FALSE, TRUE)) {
    assert_relative(
      bessel_i(-parity_x, order, scaled),
      (-1)^order * bessel_i(parity_x, order, scaled)
    )
  }
  stopifnot(
    identical(bessel_i(0, order), if (order == 0L) 1 else 0),
    identical(bessel_i(0, order, TRUE), if (order == 0L) 1 else 0),
    identical(bessel_i(Inf, order), Inf),
    identical(bessel_i(-Inf, order), (-1)^order * Inf),
    identical(bessel_i(c(-Inf, Inf), order, TRUE), c(0, 0))
  )
}
stopifnot(
  identical(bessel_i(0, 0.5), 0),
  identical(bessel_i(Inf, 0.5), Inf),
  identical(bessel_i(Inf, 0.5, TRUE), 0),
  identical(bessel_i(numeric(), 0), numeric()),
  identical(bessel_i(0:3, 1L), bessel_i(as.double(0:3), 1))
)

# A real-valued fractional-order function has no supported negative-x value.
# Report the domain problem once for the entire vector and retain valid data.
domain_warnings <- character()
domain_result <- withCallingHandlers(
  bessel_i(c(-1, -2, -Inf, 0, 1, NA_real_, NaN), 0.5),
  warning = function(w) {
    domain_warnings <<- c(domain_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
stopifnot(
  length(domain_warnings) == 1L,
  all(is.nan(domain_result[1:3])),
  identical(domain_result[4L], 0),
  is.na(domain_result[6L]),
  !is.nan(domain_result[6L]),
  is.nan(domain_result[7L])
)
assert_relative(domain_result[5L], sqrt(2 / pi) * sinh(1))

missing_result <- bessel_i(c(NA_real_, NaN, 0), 0)
stopifnot(
  is.na(missing_result[1L]),
  !is.nan(missing_result[1L]),
  is.nan(missing_result[2L]),
  identical(missing_result[3L], 1)
)

# This order exceeds the C integer range. Its value at x=1 underflows to
# zero; accepting it must not truncate the order or require order-sized memory.
stopifnot(identical(
  suppressWarnings(bessel_i(1, as.double(.Machine$integer.max) + 1, TRUE)),
  0
))

# These are public argument-contract checks, including empty x: validation
# must not disappear merely because there are no values to evaluate.
for (bad_order in list(numeric(), c(0, 1), -1, NA_real_, NaN, Inf,
                       "0", TRUE, 1i, list(0))) {
  assert_error(bessel_i(1, bad_order))
  assert_error(bessel_i(numeric(), bad_order))
}
for (bad_flag in list(logical(), c(TRUE, FALSE), NA, 0, 1, "TRUE", list(TRUE))) {
  assert_error(bessel_i(1, 0, bad_flag))
}
for (bad_x in list("1", TRUE, 1i, list(1), new.env(parent = emptyenv()))) {
  assert_error(bessel_i(bad_x, 0))
}

# The registered native interface must defend its own boundary; an R wrapper
# is not protection against callers invoking .Call directly.
native_bessel <- function(x, nu, scaled) {
  .Call("C_bessel_i", x, nu, scaled, PACKAGE = "INLAcircular")
}
stopifnot(identical(native_bessel(0:1, 0L, FALSE), bessel_i(0:1, 0)))
assert_error(native_bessel(list(1), 0, FALSE))
assert_error(native_bessel(1, numeric(), FALSE))
assert_error(native_bessel(1, c(0, 1), FALSE))
assert_error(native_bessel(1, NA_real_, FALSE))
assert_error(native_bessel(1, 0, logical()))
assert_error(native_bessel(1, 0, NA))
assert_error(native_bessel(1, 0, 1L))

# Huge orders near the exponential transition lose precision even though
# the result is representable. Return the approximation with one warning.
precision_warnings <- character()
precision_result <- withCallingHandlers(
  bessel_i(rep(6627434193491816, 2), 1e16),
  warning = function(w) {
    precision_warnings <<- c(precision_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
stopifnot(
  length(precision_warnings) == 1L,
  grepl("accuracy", precision_warnings, fixed = TRUE),
  all(is.finite(precision_result)),
  all(precision_result > 0)
)
