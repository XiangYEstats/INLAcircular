#ifndef INLACIRCULAR_BESSEL_I_H
#define INLACIRCULAR_BESSEL_I_H

/* A standalone C interface: no R or external numerical-library dependency. */
enum {
    BESSEL_I_SUCCESS = 0,
    BESSEL_I_DOMAIN = 1,
    BESSEL_I_UNDERFLOW = 2,
    BESSEL_I_OVERFLOW = 3,
    BESSEL_I_CONVERGENCE = 4,
    BESSEL_I_PRECISION = 5
};

/* Evaluate I_nu(x), or exp(-x) I_nu(x) when scaled is nonzero.
 * x and nu must be finite and nonnegative. On a range error, result is
 * zero (underflow) or positive infinity (overflow). Domain errors and
 * nonconvergence return NaN. Precision loss retains an approximate result.
 * The result pointer must not be NULL.
 */
int bessel_i(double x, double nu, int scaled, double *result);

#endif
