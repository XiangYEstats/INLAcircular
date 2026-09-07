#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <float.h>
#include <math.h>

#include "INLAcirc_common.h"

/*
 * R-facing penalized-complexity priors for circular distributions.
 *
 * This file is deliberately private to the R package.  None of these
 * routines is needed by, or shared with, the standalone INLA sources in
 * src-cloglike/.
 */

#define INLACIRC_LOG_2_LOCAL \
    0.693147180559945309417232121458176568
#define INLACIRC_PC_CARD0_MAX_DISTANCE \
    0.553942974899090755300746618256192919
#define INLACIRC_PC_CARD_MAX_DISTANCE \
    0.8325546111576977

typedef double (*INLAcirc_pc_scalar_function)(double, double);
typedef double (*INLAcirc_pc_quantile_function)(double, double, int);

static double INLAcirc_log_one_minus_exp_neg(double x)
{
    if (x == 0.0) {
        return -INFINITY;
    }
    return log(-expm1(-x));
}

/* log(1-exp(x)) for x <= 0. */
static double INLAcirc_log_one_minus_exp(double x)
{
    if (x > -INLACIRC_LOG_2_LOCAL) {
        return log(-expm1(x));
    }
    return log1p(-exp(x));
}

static double INLAcirc_cardioid_sqrt_term(double kappa)
{
    return sqrt((1.0 - 2.0 * kappa) * (1.0 + 2.0 * kappa));
}

/* d_0(kappa)^2 for the cardioid prior with uniform base. */
static double INLAcirc_pc_card0_distance_squared(double kappa)
{
    const double s = INLAcirc_cardioid_sqrt_term(kappa);
    const double t = 4.0 * kappa * kappa / (1.0 + s);

    return t + log1p(-0.5 * t);
}

static void INLAcirc_pc_card0_geometry(double kappa,
                                       double *distance,
                                       double *log_abs_derivative)
{
    if (kappa == 0.0) {
        *distance = 0.0;
        *log_abs_derivative = 0.0;
        return;
    }

    if (kappa < 1e-4) {
        const double kappa2 = kappa * kappa;
        const double distance_factor =
            1.0 + kappa2 * (0.25 + kappa2 *
                (29.0 / 96.0 + kappa2 * (211.0 / 384.0)));
        const double derivative =
            1.0 + kappa2 * (0.75 + kappa2 *
                (145.0 / 96.0 + kappa2 * (1477.0 / 384.0)));
        *distance = kappa * distance_factor;
        *log_abs_derivative = log(derivative);
        return;
    }

    {
        const double s = INLAcirc_cardioid_sqrt_term(kappa);
        const double distance_squared =
            INLAcirc_pc_card0_distance_squared(kappa);

        *distance = sqrt(fmax(0.0, distance_squared));
        *log_abs_derivative =
            log(2.0 * kappa) - log1p(s) - log(*distance);
    }
}

/* d_{1/2}(kappa)^2 near kappa=1/2, expressed as a series in s. */
static double INLAcirc_pc_card_endpoint_distance_squared(double s)
{
    const double polynomial =
        1.0 / 3.0 + s * (-1.0 / 8.0 + s *
        (1.0 / 5.0 + s * (-5.0 / 48.0 + s *
        (1.0 / 7.0 + s * (-11.0 / 128.0 + s *
        (1.0 / 9.0 + s * (-93.0 / 1280.0 + s *
        (1.0 / 11.0 + s * (-193.0 / 3072.0 + s *
        (1.0 / 13.0))))))))));

    return s * s * s * polynomial;
}

/* d_{1/2}(kappa)^2 for the cardioid prior with kappa=1/2 base. */
static double INLAcirc_pc_card_distance_squared(double kappa)
{
    const double s = INLAcirc_cardioid_sqrt_term(kappa);

    if (s < 1e-3) {
        return INLAcirc_pc_card_endpoint_distance_squared(s);
    }
    return (1.0 - 2.0 * kappa) - s + log1p(s);
}

static void INLAcirc_pc_card_geometry(double kappa,
                                      double *distance,
                                      double *log_abs_derivative)
{
    const double s = INLAcirc_cardioid_sqrt_term(kappa);

    if (kappa == 0.0) {
        *distance = INLACIRC_PC_CARD_MAX_DISTANCE;
        *log_abs_derivative = -log(*distance);
        return;
    }
    if (kappa == 0.5) {
        *distance = 0.0;
        *log_abs_derivative = INFINITY;
        return;
    }

    *distance = sqrt(fmax(0.0,
                          INLAcirc_pc_card_distance_squared(kappa)));
    *log_abs_derivative =
        log((1.0 - 2.0 * kappa) + s) - log1p(s) - log(*distance);
}

/* Stable value of d(0)-d(kappa), used to invert close to the boundary atom. */
static double INLAcirc_pc_card_distance_deficit(double kappa)
{
    double distance;
    double unused_log_abs_derivative;
    double difference_of_squares;

    if (kappa == 0.0) {
        return 0.0;
    }
    if (kappa == 0.5) {
        return INLACIRC_PC_CARD_MAX_DISTANCE;
    }

    INLAcirc_pc_card_geometry(kappa, &distance,
                              &unused_log_abs_derivative);
    difference_of_squares =
        2.0 * kappa - INLAcirc_pc_card0_distance_squared(kappa);
    return difference_of_squares /
           (INLACIRC_PC_CARD_MAX_DISTANCE + distance);
}

static void INLAcirc_pc_wc_geometry(double kappa,
                                    double *distance,
                                    double *log_abs_derivative)
{
    if (kappa == 0.0) {
        *distance = 0.0;
        *log_abs_derivative = 0.0;
        return;
    }
    if (kappa == 1.0) {
        *distance = INFINITY;
        *log_abs_derivative = INFINITY;
        return;
    }

    if (kappa < 1e-4) {
        const double kappa2 = kappa * kappa;
        const double distance_factor =
            1.0 + kappa2 * (0.25 + kappa2 *
                (13.0 / 96.0 + kappa2 * (35.0 / 384.0)));
        const double derivative =
            1.0 + kappa2 * (0.75 + kappa2 *
                (65.0 / 96.0 + kappa2 * (245.0 / 384.0)));

        *distance = kappa * distance_factor;
        *log_abs_derivative = log(derivative);
        return;
    }

    {
        const double log_one_minus_kappa_squared =
            log1p(-kappa) + log1p(kappa);

        *distance = sqrt(-log_one_minus_kappa_squared);
        *log_abs_derivative =
            log(kappa) - log_one_minus_kappa_squared - log(*distance);
    }
}

static double INLAcirc_pc_card0_log_density(double kappa, double lambda)
{
    double distance;
    double log_abs_derivative;
    double lambda_max_distance;
    double log_normalizer;

    if (isnan(kappa) || isnan(lambda) || !isfinite(lambda) ||
        lambda <= 0.0) {
        return NAN;
    }
    if (kappa < 0.0 || kappa > 0.5) {
        return -INFINITY;
    }

    INLAcirc_pc_card0_geometry(kappa, &distance, &log_abs_derivative);
    lambda_max_distance = lambda * INLACIRC_PC_CARD0_MAX_DISTANCE;
    if (lambda_max_distance < DBL_MIN) {
        return -log(INLACIRC_PC_CARD0_MAX_DISTANCE) +
               log_abs_derivative;
    }
    log_normalizer = INLAcirc_log_one_minus_exp_neg(lambda_max_distance);

    return log(lambda) - lambda * distance - log_normalizer +
           log_abs_derivative;
}

static double INLAcirc_pc_card0_log_cdf(double kappa, double lambda)
{
    double distance;
    double unused_log_abs_derivative;
    double lambda_max_distance;
    double log_normalizer;

    if (isnan(kappa) || isnan(lambda) || !isfinite(lambda) ||
        lambda <= 0.0) {
        return NAN;
    }
    if (kappa < 0.0) {
        return -INFINITY;
    }
    if (kappa >= 0.5) {
        return 0.0;
    }
    if (kappa == 0.0) {
        return -INFINITY;
    }

    INLAcirc_pc_card0_geometry(kappa, &distance,
                               &unused_log_abs_derivative);
    lambda_max_distance = lambda * INLACIRC_PC_CARD0_MAX_DISTANCE;
    if (lambda_max_distance < DBL_MIN) {
        return log(distance) - log(INLACIRC_PC_CARD0_MAX_DISTANCE);
    }
    log_normalizer = INLAcirc_log_one_minus_exp_neg(lambda_max_distance);

    {
        const double lambda_distance = lambda * distance;
        const double log_numerator = (lambda_distance < DBL_MIN)
            ? log(lambda) + log(distance)
            : INLAcirc_log_one_minus_exp_neg(lambda_distance);
        return log_numerator - log_normalizer;
    }
}

static double INLAcirc_pc_card0_quantile(double supplied_probability,
                                         double lambda,
                                         int probability_is_log)
{
    const double lower_log_kappa = log(nextafter(0.0, 1.0));
    const double upper_log_kappa = log(0.5);
    const double lambda_max_distance =
        lambda * INLACIRC_PC_CARD0_MAX_DISTANCE;
    double probability;
    double log_normalizer;
    double target_distance;
    double target_log_distance;
    double lower = lower_log_kappa;
    double upper = upper_log_kappa;

    if (isnan(supplied_probability) || isnan(lambda) ||
        !isfinite(lambda) || lambda <= 0.0 ||
        (probability_is_log && supplied_probability > 0.0) ||
        (!probability_is_log &&
         (supplied_probability < 0.0 || supplied_probability > 1.0))) {
        return NAN;
    }

    if (probability_is_log) {
        if (supplied_probability == -INFINITY) {
            return 0.0;
        }
        if (supplied_probability == 0.0) {
            return 0.5;
        }
        if (lambda_max_distance < DBL_MIN) {
            target_distance = exp(
                supplied_probability +
                log(INLACIRC_PC_CARD0_MAX_DISTANCE));
        } else {
            log_normalizer =
                INLAcirc_log_one_minus_exp_neg(lambda_max_distance);
            target_distance = -INLAcirc_log_one_minus_exp(
                supplied_probability + log_normalizer) / lambda;
        }
    } else {
        probability = supplied_probability;
        if (probability == 0.0) {
            return 0.0;
        }
        if (probability == 1.0) {
            return 0.5;
        }
        if (lambda_max_distance < DBL_MIN) {
            target_distance =
                probability * INLACIRC_PC_CARD0_MAX_DISTANCE;
        } else {
            const double normalizer = -expm1(-lambda_max_distance);
            target_distance = -log1p(-probability * normalizer) / lambda;
        }
    }

    if (target_distance <= 0.0) {
        return 0.0;
    }
    if (target_distance >= INLACIRC_PC_CARD0_MAX_DISTANCE) {
        return 0.5;
    }
    target_log_distance = log(target_distance);

    for (int iteration = 0; iteration < 120; ++iteration) {
        const double middle = lower + 0.5 * (upper - lower);
        const double kappa = exp(middle);
        double distance;
        double unused_log_abs_derivative;

        INLAcirc_pc_card0_geometry(kappa, &distance,
                                   &unused_log_abs_derivative);
        if (log(distance) < target_log_distance) {
            lower = middle;
        } else {
            upper = middle;
        }
    }

    return exp(lower + 0.5 * (upper - lower));
}

static double INLAcirc_pc_card_log_density(double kappa, double lambda)
{
    double distance;
    double log_abs_derivative;

    if (isnan(kappa) || isnan(lambda) || !isfinite(lambda) ||
        lambda <= 0.0) {
        return NAN;
    }
    if (kappa < 0.0 || kappa > 0.5) {
        return -INFINITY;
    }
    if (kappa == 0.5) {
        return INFINITY;
    }

    INLAcirc_pc_card_geometry(kappa, &distance, &log_abs_derivative);
    return log(lambda) - lambda * distance + log_abs_derivative;
}

static double INLAcirc_pc_card_log_cdf(double kappa, double lambda)
{
    double distance;
    double unused_log_abs_derivative;

    if (isnan(kappa) || isnan(lambda) || !isfinite(lambda) ||
        lambda <= 0.0) {
        return NAN;
    }
    if (kappa < 0.0) {
        return -INFINITY;
    }
    if (kappa >= 0.5) {
        return 0.0;
    }

    INLAcirc_pc_card_geometry(kappa, &distance,
                              &unused_log_abs_derivative);
    return -lambda * distance;
}

static double INLAcirc_pc_card_quantile(double supplied_probability,
                                        double lambda,
                                        int probability_is_log)
{
    const double atom_log_probability =
        -lambda * INLACIRC_PC_CARD_MAX_DISTANCE;
    const double lower_log_kappa = log(nextafter(0.0, 1.0));
    const double upper_log_kappa = log(0.5);
    double log_probability;
    double target_deficit;
    double lower = lower_log_kappa;
    double upper = upper_log_kappa;

    if (isnan(supplied_probability) || isnan(lambda) ||
        !isfinite(lambda) || lambda <= 0.0 ||
        (probability_is_log && supplied_probability > 0.0) ||
        (!probability_is_log &&
         (supplied_probability < 0.0 || supplied_probability > 1.0))) {
        return NAN;
    }

    if (probability_is_log) {
        log_probability = supplied_probability;
        if (log_probability <= atom_log_probability) {
            return 0.0;
        }
    } else {
        const double atom_probability = exp(atom_log_probability);
        if (supplied_probability <= atom_probability) {
            return 0.0;
        }
        log_probability = log(supplied_probability);
    }
    if (log_probability == 0.0) {
        return 0.5;
    }

    target_deficit =
        (log_probability - atom_log_probability) / lambda;
    if (target_deficit <= 0.0) {
        return 0.0;
    }
    if (target_deficit >= INLACIRC_PC_CARD_MAX_DISTANCE) {
        return 0.5;
    }

    for (int iteration = 0; iteration < 120; ++iteration) {
        const double middle = lower + 0.5 * (upper - lower);
        const double kappa = exp(middle);
        const double deficit =
            INLAcirc_pc_card_distance_deficit(kappa);

        if (deficit < target_deficit) {
            lower = middle;
        } else {
            upper = middle;
        }
    }

    return exp(lower + 0.5 * (upper - lower));
}

static double INLAcirc_pc_wc_log_density(double kappa, double lambda)
{
    double distance;
    double log_abs_derivative;

    if (isnan(kappa) || isnan(lambda) || !isfinite(lambda) ||
        lambda <= 0.0) {
        return NAN;
    }
    if (kappa < 0.0 || kappa > 1.0) {
        return -INFINITY;
    }
    if (kappa == 1.0) {
        return INFINITY;
    }

    INLAcirc_pc_wc_geometry(kappa, &distance, &log_abs_derivative);
    return log(lambda) - lambda * distance + log_abs_derivative;
}

static double INLAcirc_pc_wc_log_cdf(double kappa, double lambda)
{
    double distance;
    double unused_log_abs_derivative;

    if (isnan(kappa) || isnan(lambda) || !isfinite(lambda) ||
        lambda <= 0.0) {
        return NAN;
    }
    if (kappa < 0.0) {
        return -INFINITY;
    }
    if (kappa >= 1.0) {
        return 0.0;
    }
    if (kappa == 0.0) {
        return -INFINITY;
    }

    INLAcirc_pc_wc_geometry(kappa, &distance,
                            &unused_log_abs_derivative);
    {
        const double lambda_distance = lambda * distance;
        return (lambda_distance < DBL_MIN)
            ? log(lambda) + log(distance)
            : INLAcirc_log_one_minus_exp_neg(lambda_distance);
    }
}

static double INLAcirc_pc_wc_quantile(double supplied_probability,
                                      double lambda,
                                      int probability_is_log)
{
    double log_survival;
    double distance;
    double distance_squared;

    if (isnan(supplied_probability) || isnan(lambda) ||
        !isfinite(lambda) || lambda <= 0.0 ||
        (probability_is_log && supplied_probability > 0.0) ||
        (!probability_is_log &&
         (supplied_probability < 0.0 || supplied_probability > 1.0))) {
        return NAN;
    }

    if (probability_is_log) {
        if (supplied_probability == -INFINITY) {
            return 0.0;
        }
        if (supplied_probability == 0.0) {
            return 1.0;
        }
        log_survival = INLAcirc_log_one_minus_exp(supplied_probability);
    } else {
        if (supplied_probability == 0.0) {
            return 0.0;
        }
        if (supplied_probability == 1.0) {
            return 1.0;
        }
        log_survival = log1p(-supplied_probability);
    }

    distance = -log_survival / lambda;
    if (distance < sqrt(DBL_MIN)) {
        /*
         * For this range, distance^2 underflows (or is subnormal), while
         * sqrt(1 - exp(-distance^2)) rounds to distance itself.
         */
        return distance;
    }
    distance_squared = distance * distance;
    if (!isfinite(distance_squared)) {
        return nextafter(1.0, 0.0);
    }
    {
        const double quantile = sqrt(-expm1(-distance_squared));
        return (quantile < 1.0) ? quantile : nextafter(1.0, 0.0);
    }
}

static SEXP INLAcirc_pc_apply(SEXP kappa,
                              SEXP lambda,
                              SEXP log_probability,
                              INLAcirc_pc_scalar_function function)
{
    const int count = Rf_length(kappa);
    const double *concentrations = REAL(kappa);
    const double *rates = REAL(lambda);
    const int return_log = Rf_asLogical(log_probability);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        const double log_value = function(concentrations[i], rates[i]);
        result[i] = return_log ? log_value : exp(log_value);
    }

    UNPROTECT(1);
    return output;
}

static SEXP INLAcirc_pc_quantile_apply(
    SEXP p,
    SEXP lambda,
    SEXP log_probability,
    INLAcirc_pc_quantile_function function)
{
    const int count = Rf_length(p);
    const double *probabilities = REAL(p);
    const double *rates = REAL(lambda);
    const int probability_is_log = Rf_asLogical(log_probability);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        result[i] = function(probabilities[i], rates[i],
                             probability_is_log);
    }

    UNPROTECT(1);
    return output;
}

static SEXP INLAcirc_pc_random_apply(
    SEXP n,
    SEXP lambda,
    INLAcirc_pc_quantile_function function)
{
    const int count = Rf_asInteger(n);
    const double *rates = REAL(lambda);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    GetRNGstate();
    for (int i = 0; i < count; ++i) {
        result[i] = function(unif_rand(), rates[i], 0);
    }
    PutRNGstate();

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_dpc_vm0(SEXP kappa,
                        SEXP lambda,
                        SEXP log_probability)
{
    const int count = Rf_length(kappa);
    const double *concentrations = REAL(kappa);
    const double *rates = REAL(lambda);
    const int return_log = Rf_asLogical(log_probability);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        const double log_density =
            INLAcirc_pc_vm0_log_density(concentrations[i], rates[i]);
        result[i] = return_log ? log_density : exp(log_density);
    }

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_ppc_vm0(SEXP q,
                        SEXP lambda,
                        SEXP log_probability)
{
    const int count = Rf_length(q);
    const double *quantiles = REAL(q);
    const double *rates = REAL(lambda);
    const int return_log = Rf_asLogical(log_probability);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        const double log_cdf =
            INLAcirc_pc_vm0_log_cdf(quantiles[i], rates[i]);
        result[i] = return_log ? log_cdf : exp(log_cdf);
    }

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_qpc_vm0(SEXP p,
                        SEXP lambda,
                        SEXP log_probability)
{
    const int count = Rf_length(p);
    const double *probabilities = REAL(p);
    const double *rates = REAL(lambda);
    const int probabilities_are_log = Rf_asLogical(log_probability);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        const double supplied_probability = probabilities[i];
        double probability;

        if (isnan(supplied_probability) ||
            (probabilities_are_log && supplied_probability > 0.0) ||
            (!probabilities_are_log &&
             (supplied_probability < 0.0 || supplied_probability > 1.0))) {
            result[i] = NAN;
            continue;
        }

        if (probabilities_are_log) {
            probability = exp(supplied_probability);
        } else {
            probability = supplied_probability;
        }
        result[i] = INLAcirc_pc_vm0_quantile(probability, rates[i]);
    }

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_rpc_vm0(SEXP n, SEXP lambda)
{
    const int count = Rf_asInteger(n);
    const double *rates = REAL(lambda);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    GetRNGstate();
    for (int i = 0; i < count; ++i) {
        result[i] = INLAcirc_pc_vm0_quantile(unif_rand(), rates[i]);
    }
    PutRNGstate();

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_dpc_vminf(SEXP kappa,
                          SEXP lambda,
                          SEXP log_probability)
{
    const int count = Rf_length(kappa);
    const double *concentrations = REAL(kappa);
    const double *rates = REAL(lambda);
    const int return_log = Rf_asLogical(log_probability);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        const double log_density =
            INLAcirc_pc_vminf_log_density(concentrations[i], rates[i]);
        result[i] = return_log ? log_density : exp(log_density);
    }

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_ppc_vminf(SEXP q,
                          SEXP lambda,
                          SEXP log_probability)
{
    const int count = Rf_length(q);
    const double *quantiles = REAL(q);
    const double *rates = REAL(lambda);
    const int return_log = Rf_asLogical(log_probability);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        const double log_cdf =
            INLAcirc_pc_vminf_log_cdf(quantiles[i], rates[i]);
        result[i] = return_log ? log_cdf : exp(log_cdf);
    }

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_qpc_vminf(SEXP p,
                          SEXP lambda,
                          SEXP log_probability)
{
    const int count = Rf_length(p);
    const double *probabilities = REAL(p);
    const double *rates = REAL(lambda);
    const int probabilities_are_log = Rf_asLogical(log_probability);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        const double supplied_probability = probabilities[i];
        double probability;

        if (isnan(supplied_probability) ||
            (probabilities_are_log && supplied_probability > 0.0) ||
            (!probabilities_are_log &&
             (supplied_probability < 0.0 || supplied_probability > 1.0))) {
            result[i] = NAN;
            continue;
        }

        probability = probabilities_are_log
                          ? exp(supplied_probability)
                          : supplied_probability;
        result[i] = INLAcirc_pc_vminf_quantile(probability, rates[i]);
    }

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_rpc_vminf(SEXP n, SEXP lambda)
{
    const int count = Rf_asInteger(n);
    const double *rates = REAL(lambda);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    GetRNGstate();
    for (int i = 0; i < count; ++i) {
        result[i] = INLAcirc_pc_vminf_quantile(unif_rand(), rates[i]);
    }
    PutRNGstate();

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_dpc_card0(SEXP kappa, SEXP lambda, SEXP log_probability)
{
    return INLAcirc_pc_apply(kappa, lambda, log_probability,
                             INLAcirc_pc_card0_log_density);
}

SEXP INLAcirc_C_ppc_card0(SEXP q, SEXP lambda, SEXP log_probability)
{
    return INLAcirc_pc_apply(q, lambda, log_probability,
                             INLAcirc_pc_card0_log_cdf);
}

SEXP INLAcirc_C_qpc_card0(SEXP p, SEXP lambda, SEXP log_probability)
{
    return INLAcirc_pc_quantile_apply(p, lambda, log_probability,
                                      INLAcirc_pc_card0_quantile);
}

SEXP INLAcirc_C_rpc_card0(SEXP n, SEXP lambda)
{
    return INLAcirc_pc_random_apply(n, lambda,
                                    INLAcirc_pc_card0_quantile);
}

SEXP INLAcirc_C_dpc_card(SEXP kappa, SEXP lambda, SEXP log_probability)
{
    return INLAcirc_pc_apply(kappa, lambda, log_probability,
                             INLAcirc_pc_card_log_density);
}

SEXP INLAcirc_C_ppc_card(SEXP q, SEXP lambda, SEXP log_probability)
{
    return INLAcirc_pc_apply(q, lambda, log_probability,
                             INLAcirc_pc_card_log_cdf);
}

SEXP INLAcirc_C_qpc_card(SEXP p, SEXP lambda, SEXP log_probability)
{
    return INLAcirc_pc_quantile_apply(p, lambda, log_probability,
                                      INLAcirc_pc_card_quantile);
}

SEXP INLAcirc_C_rpc_card(SEXP n, SEXP lambda)
{
    return INLAcirc_pc_random_apply(n, lambda,
                                    INLAcirc_pc_card_quantile);
}

SEXP INLAcirc_C_dpc_wc(SEXP kappa, SEXP lambda, SEXP log_probability)
{
    return INLAcirc_pc_apply(kappa, lambda, log_probability,
                             INLAcirc_pc_wc_log_density);
}

SEXP INLAcirc_C_ppc_wc(SEXP q, SEXP lambda, SEXP log_probability)
{
    return INLAcirc_pc_apply(q, lambda, log_probability,
                             INLAcirc_pc_wc_log_cdf);
}

SEXP INLAcirc_C_qpc_wc(SEXP p, SEXP lambda, SEXP log_probability)
{
    return INLAcirc_pc_quantile_apply(p, lambda, log_probability,
                                      INLAcirc_pc_wc_quantile);
}

SEXP INLAcirc_C_rpc_wc(SEXP n, SEXP lambda)
{
    return INLAcirc_pc_random_apply(n, lambda,
                                    INLAcirc_pc_wc_quantile);
}
