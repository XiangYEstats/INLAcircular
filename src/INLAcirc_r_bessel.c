#include <R.h>
#include <Rinternals.h>

#include "INLAcirc_common.h"

SEXP INLAcirc_C_bessel_i(SEXP x, SEXP nu, SEXP expon_scaled)
{
    const int length_x = Rf_length(x);
    const double *values = REAL(x);
    const double order = Rf_asReal(nu);
    const int scaled = Rf_asLogical(expon_scaled);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, length_x));
    double *result = REAL(output);

    for (int i = 0; i < length_x; ++i) {
        result[i] = ISNAN(values[i])
                        ? NA_REAL
                        : INLAcirc_bessel_i(values[i], order, scaled);
    }

    UNPROTECT(1);
    return output;
}
