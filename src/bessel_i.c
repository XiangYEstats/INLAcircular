/*
 * Standalone modified Bessel I, adapted from GNU Scientific Library 2.8:
 *   specfunc/bessel_Inu.c, Copyright (C) 1996-2000 Gerard Jungman
 *   specfunc/bessel.c,     Copyright (C) 1996-2003 Gerard Jungman
 * Original author: G. Jungman.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <https://www.gnu.org/licenses/>.
 *
 * Adaptations for INLAcircular:
 * - retain the GSL Taylor / uniform-asymptotic / continued-fraction
 *   branch strategy for finite nonnegative arguments and orders;
 * - replace GSL gamma helpers with a logarithmic C99 lgamma/Stirling
 *   prefactor and evaluate the final exponential only once;
 * - rearrange the Debye polynomials in inverse hypot(x,nu), avoiding
 *   divisions by nu and allowing order zero in the uniform expansion;
 * - omit the unreachable x < 2 Temme branch (the Taylor condition is
 *   always satisfied there for nu >= 0);
 * - replace GSL errors and result structures with explicit status codes;
 * - reject nonconvergence and avoid intermediate overflow/underflow.
 *
 * No R API, GSL library, allocation, or mutable global state is used.
 */

#include "bessel_i.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

static const long double bessel_i_log_two =
    0.693147180559945309417232121458L;
static const long double bessel_i_log_sqrt_two_pi =
    0.918938533204672741780329736406L;

/* (x/2)^nu / Gamma(nu+1), in logarithmic form. The Stirling form
 * avoids subtracting two large quantities, or overflowing lgamma,
 * when the order is large. At nu >= 32 the omitted term is < 1e-21.
 */
static long double bessel_i_log_prefactor(double x, double nu)
{
    const long double log_half_x = logl(x) - bessel_i_log_two;

    if (nu == 0.0) {
        return 0.0;
    }
    if (nu < 32.0) {
        return nu * log_half_x - lgammal((long double) nu + 1.0L);
    }
    {
        const double inverse = 1.0 / nu;
        const double inverse2 = inverse * inverse;
        const double correction = inverse *
            (1.0 / 12.0 + inverse2 *
             (-1.0 / 360.0 + inverse2 *
              (1.0 / 1260.0 + inverse2 *
               (-1.0 / 1680.0 + inverse2 *
                (1.0 / 1188.0 - inverse2 * 691.0 / 360360.0)))));
        const long double log_order = logl(nu);
        return nu * (log_half_x - log_order + 1.0) -
            0.5 * log_order - bessel_i_log_sqrt_two_pi - correction;
    }
}

/* Positive Taylor series from gsl_sf_bessel_IJ_taylor_e, sign = +1.
 * This branch has x^2 / (4 (nu+1)) < 2.5, so its terms remain bounded.
 */
static int bessel_i_taylor_log(double x, double nu, int scaled,
                             long double *log_result)
{
    double term = 1.0;
    double sum = 1.0;
    int k;

    for (k = 1; k <= 100; ++k) {
        term *= 0.25 * (x / (nu + k)) * (x / k);
        sum += term;
        if (term < DBL_EPSILON * sum) {
            *log_result = bessel_i_log_prefactor(x, nu) + logl(sum) -
                (scaled ? x : 0.0);
            return BESSEL_I_SUCCESS;
        }
    }
    return BESSEL_I_CONVERGENCE;
}

/* GSL gsl_sf_bessel_I_CF1_ser: I_(nu+1)(x) / I_nu(x), using
 * Gautschi's (Euler) equivalent series for the continued fraction.
 */
static int bessel_i_ratio(double x, double nu, double *ratio)
{
    double term = 1.0;
    double sum = 1.0;
    double rho = 0.0;
    int k;

    for (k = 1; k < 20000; ++k) {
        const double a = 0.25 * (x / (nu + k)) * x / (nu + k + 1.0);
        rho = -a * (1.0 + rho) / (1.0 + a * (1.0 + rho));
        term *= rho;
        sum += term;
        if (fabs(term) < DBL_EPSILON * fabs(sum)) {
            *ratio = x / (2.0 * (nu + 1.0)) * sum;
            return BESSEL_I_SUCCESS;
        }
    }
    return BESSEL_I_CONVERGENCE;
}

/* GSL gsl_sf_bessel_K_scaled_steed_temme_CF2: Thompson-Barnett-
 * Temme CF2 for exp(x) K_mu(x) and exp(x) K_(mu+1)(x).
 * The caller guarantees x >= 2 and |mu| <= 1/2.
 */
static int bessel_i_scaled_k_pair(double x, double mu,
                                  double *k_mu, double *k_mu_plus_one)
{
    double b = 2.0 * (1.0 + x);
    double d = 1.0 / b;
    double delta_h = d;
    double h = d;
    double q = 0.0;
    double q_next = 1.0;
    double a = -(0.25 - mu * mu);
    const double a_first = a;
    double c = -a;
    double big_q = -a;
    double sum = 1.0 + big_q * delta_h;
    int i;

    for (i = 2; i <= 10000; ++i) {
        double temporary;
        double delta_sum;
        a -= 2.0 * (i - 1);
        c = -a * c / i;
        temporary = (q - b * q_next) / a;
        q = q_next;
        q_next = temporary;
        big_q += c * q_next;
        b += 2.0;
        d = 1.0 / (b + a * d);
        delta_h = (b * d - 1.0) * delta_h;
        h += delta_h;
        delta_sum = big_q * delta_h;
        sum += delta_sum;
        if (fabs(delta_sum) < DBL_EPSILON * fabs(sum)) {
            const double half_pi = 1.570796326794896619231321691640;
            h *= -a_first;
            *k_mu = sqrt(half_pi / x) / sum;
            *k_mu_plus_one = *k_mu * (mu + x + 0.5 - h) / x;
            return BESSEL_I_SUCCESS;
        }
    }
    return BESSEL_I_CONVERGENCE;
}

/* GSL's five-term uniform Debye expansion (A&S 9.7.7), rewritten
 * algebraically with r = hypot(x,nu), t = nu/r and u_k(t)/nu^k.
 * This form remains defined at nu = 0 and when x/nu would overflow.
 */
static long double bessel_i_uniform_log(double x, double nu, int scaled,
                                       int *precision_loss)
{
    const double scale = fmax(x, nu);
    const double small_ratio = fmin(x, nu) / scale;
    const double root = hypot(1.0, small_ratio);
    const double inverse_r = (1.0 / scale) / root;
    const double t = (nu / scale) / root;
    const double t2 = t * t;
    const long double log_r = logl(scale) + logl(root);
    const double u1 = (3.0 - 5.0 * t2) / 24.0;
    const double u2 = (81.0 + t2 * (-462.0 + 385.0 * t2)) / 1152.0;
    const double u3 = (30375.0 + t2 * (-369603.0 + t2 *
                      (765765.0 - 425425.0 * t2))) / 414720.0;
    const double u4 = (4465125.0 + t2 * (-94121676.0 + t2 *
                      (349922430.0 + t2 * (-446185740.0 +
                       185910725.0 * t2)))) / 39813120.0;
    const double u5 = (1519035525.0 + t2 * (-49286948607.0 + t2 *
                      (284499769554.0 + t2 * (-614135872350.0 + t2 *
                       (566098157625.0 - 188699385875.0 * t2))))) /
                      6688604160.0;
    const double correction = inverse_r * (u1 + inverse_r *
        (u2 + inverse_r * (u3 + inverse_r * (u4 + inverse_r * u5))));
    long double exponent;

    if (nu <= x) {
        const long double ratio = (long double) nu / x;
        const long double local_root = hypotl(1.0L, ratio);
        /* r - x = nu * ratio/(sqrt(1+ratio^2)+1). */
        const long double scaled_exponent =
            nu * (ratio / (local_root + 1.0L) - asinhl(ratio));
        exponent = scaled ? scaled_exponent : x + scaled_exponent;
    } else {
        const long double ratio = (long double) x / nu;
        const long double local_root = hypotl(1.0L, ratio);
        const long double logarithm = log1pl(local_root) - logl(ratio);
        exponent = nu * (scaled ?
            1.0L / (local_root + ratio) - logarithm :
            local_root - logarithm);
        if (!scaled) {
            /* Around x/nu = 0.662743..., the unscaled exponent is the
             * difference of nearly equal terms multiplied by nu.
             * Long-double arithmetic reduces, but cannot eliminate,
             * the loss for very large orders. Estimate its absolute
             * log error, allowing for the elementary operations above.
             * This is a roundoff estimate, not a rigorous error bound.
             */
            const long double log_error = 8.0L * LDBL_EPSILON * nu *
                (fabsl(local_root) + fabsl(logarithm));
            *precision_loss = log_error > 1.0e-10L;
        }
    }

    return exponent - bessel_i_log_sqrt_two_pi - 0.5 * log_r +
        log1pl(correction);
}

/* GSL's intermediate-argument branch, normalizing I using its ratio
 * and the Wronskian with K. The preceding branches bound nu below 288,
 * so the integer conversion and forward recurrence are both bounded.
 */
static int bessel_i_cf_log(double x, double nu, int scaled,
                          long double *log_result)
{
    const int n = (int) (nu + 0.5);
    const double mu = nu - n;
    double k_current;
    double k_next;
    double ratio;
    int j;
    int status = bessel_i_scaled_k_pair(x, mu, &k_current, &k_next);

    if (status != BESSEL_I_SUCCESS) {
        return status;
    }
    for (j = 0; j < n; ++j) {
        const double previous = k_current;
        k_current = k_next;
        k_next = 2.0 * (mu + j + 1.0) / x * k_current + previous;
    }
    status = bessel_i_ratio(x, nu, &ratio);
    if (status != BESSEL_I_SUCCESS) {
        return status;
    }
    *log_result = -logl(x) - logl(k_next + ratio * k_current) +
        (scaled ? 0.0 : x);
    return BESSEL_I_SUCCESS;
}

int bessel_i(double x, double nu, int scaled, double *result)
{
    long double log_result;
    int status;
    int precision_loss = 0;

    if (result == NULL) {
        return BESSEL_I_DOMAIN;
    }
    *result = NAN;
    if (!isfinite(x) || !isfinite(nu) || x < 0.0 || nu < 0.0) {
        return BESSEL_I_DOMAIN;
    }
    if (x == 0.0) {
        *result = nu == 0.0 ? 1.0 : 0.0;
        return BESSEL_I_SUCCESS;
    }

    if (x / sqrt(nu + 1.0) < 3.162277660168379331998893544433) {
        status = bessel_i_taylor_log(x, nu, scaled, &log_result);
    } else if (hypot(x, nu) > sqrt(0.5 / cbrt(DBL_EPSILON))) {
        log_result = bessel_i_uniform_log(x, nu, scaled, &precision_loss);
        status = BESSEL_I_SUCCESS;
    } else {
        status = bessel_i_cf_log(x, nu, scaled, &log_result);
    }
    if (status != BESSEL_I_SUCCESS || isnan(log_result)) {
        return BESSEL_I_CONVERGENCE;
    }
    if (log_result > logl(DBL_MAX)) {
        *result = INFINITY;
        return BESSEL_I_OVERFLOW;
    }
    *result = (double) expl(log_result);
    if (*result == 0.0) {
        return BESSEL_I_UNDERFLOW;
    }
    if (!isfinite(*result)) {
        return BESSEL_I_OVERFLOW;
    }
    return precision_loss ? BESSEL_I_PRECISION : BESSEL_I_SUCCESS;
}
