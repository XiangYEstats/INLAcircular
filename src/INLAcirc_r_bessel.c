#include <R.h>
#include <Rinternals.h>
#include <math.h>

#include "bessel_i.h"

SEXP C_bessel_i(SEXP x, SEXP nu, SEXP expon_scaled)
{
    if ((TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP) || Rf_isFactor(x)) {
        Rf_error("'x' must be a numeric vector.");
    }
    if ((TYPEOF(nu) != REALSXP && TYPEOF(nu) != INTSXP) ||
        Rf_isFactor(nu) || XLENGTH(nu) != 1) {
        Rf_error("'nu' must be one finite, nonnegative numeric value.");
    }
    const double order = Rf_asReal(nu);
    if (!R_FINITE(order) || order < 0.0) {
        Rf_error("'nu' must be one finite, nonnegative numeric value.");
    }
    if (TYPEOF(expon_scaled) != LGLSXP || XLENGTH(expon_scaled) != 1 ||
        LOGICAL(expon_scaled)[0] == NA_LOGICAL) {
        Rf_error("'expon.scaled' must be TRUE or FALSE.");
    }

    const int scaled = LOGICAL(expon_scaled)[0];
    const int integer_order = order == floor(order);
    const double negative_sign = fmod(order, 2.0) == 0.0 ? 1.0 : -1.0;
    const R_xlen_t length_x = XLENGTH(x);
    SEXP numeric_x = PROTECT(Rf_coerceVector(x, REALSXP));
    const double *values = REAL(numeric_x);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, length_x));
    double *result = REAL(output);
    int domain_warning = 0;
    int convergence_warning = 0;
    int precision_warning = 0;

    for (R_xlen_t i = 0; i < length_x; ++i) {
        if (i % 16384 == 0) {
            R_CheckUserInterrupt();
        }
        const double value = values[i];
        if (ISNAN(value)) {
            result[i] = value;
            continue;
        }
        if (value < 0.0 && !integer_order) {
            result[i] = R_NaN;
            domain_warning = 1;
            continue;
        }
        const double sign = value < 0.0 ? negative_sign : 1.0;
        if (!R_FINITE(value)) {
            result[i] = sign * (scaled ? 0.0 : R_PosInf);
            continue;
        }
        const int status = bessel_i(fabs(value), order, scaled, &result[i]);
        result[i] *= sign;
        if (status == BESSEL_I_DOMAIN) {
            domain_warning = 1;
        } else if (status == BESSEL_I_CONVERGENCE) {
            convergence_warning = 1;
        } else if (status == BESSEL_I_PRECISION) {
            precision_warning = 1;
        }
    }

    if (domain_warning) {
        Rf_warning("NaNs produced: negative 'x' requires an integer 'nu'.");
    }
    if (convergence_warning) {
        Rf_warning("Bessel I calculation did not converge for some inputs.");
    }
    if (precision_warning) {
        Rf_warning("Bessel I accuracy may be reduced for very large orders.");
    }
    UNPROTECT(2);
    return output;
}
