#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <math.h>

#include "../src-cloglike/INLAcirc_common.h"

static double INLAcirc_approx_linear(double x,
                                     const double *grid_x,
                                     const double *grid_y,
                                     int grid_length)
{
    int left;
    int right;

    if (x <= grid_x[0]) {
        return grid_y[0];
    }
    if (x >= grid_x[grid_length - 1]) {
        return grid_y[grid_length - 1];
    }

    left = 0;
    right = grid_length - 1;
    while (left <= right) {
        const int middle = left + (right - left) / 2;

        if (grid_x[middle] == x) {
            return grid_y[middle];
        }
        if (grid_x[middle] < x) {
            left = middle + 1;
        } else {
            right = middle - 1;
        }
    }

    {
        const double x0 = grid_x[right];
        const double x1 = grid_x[left];
        const double y0 = grid_y[right];
        const double y1 = grid_y[left];
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
    }
}

SEXP INLAcirc_C_dvm(SEXP x, SEXP mu, SEXP kappa, SEXP log_probability)
{
    const int count = Rf_length(x);
    const double *observations = REAL(x);
    const double *locations = REAL(mu);
    const double *concentrations = REAL(kappa);
    const int return_log = Rf_asLogical(log_probability);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        const double concentration = concentrations[i];
        const double sine = sin((observations[i] - locations[i]) / 2.0);
        const double log_density =
            -2.0 * concentration * sine * sine - INLACIRC_LOG_2PI -
            INLAcirc_log_bessel_i0_scaled(concentration);

        result[i] = return_log ? log_density : exp(log_density);
    }

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_dlavm(SEXP x,
                      SEXP eta,
                      SEXP kappa,
                      SEXP log_probability)
{
    const int count = Rf_length(x);
    const double *observations = REAL(x);
    const double *predictors = REAL(eta);
    const double *concentrations = REAL(kappa);
    const int return_log = Rf_asLogical(log_probability);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        const double tangent_half_x = tan(observations[i] / 2.0);
        const double shifted = tangent_half_x - predictors[i];
        const double transformed = 2.0 * atan(shifted);
        const double sine = sin(transformed / 2.0);
        const double log_jacobian =
            log1p(tangent_half_x * tangent_half_x) -
            log1p(shifted * shifted);
        const double log_density =
            -2.0 * concentrations[i] * sine * sine - INLACIRC_LOG_2PI -
            INLAcirc_log_bessel_i0_scaled(concentrations[i]) +
            log_jacobian;

        result[i] = return_log ? log_density : exp(log_density);
    }

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

SEXP INLAcirc_C_pvm(SEXP q,
                    SEXP mu,
                    SEXP strategy,
                    SEXP log_probability,
                    SEXP grid_x,
                    SEXP grid_y)
{
    const int count = Rf_length(q);
    const int grid_length = Rf_length(grid_x);
    const double *quantiles = REAL(q);
    const double *locations = REAL(mu);
    const double *x_grid = REAL(grid_x);
    const double *y_grid = REAL(grid_y);
    const int strategy_code = Rf_asInteger(strategy);
    const int return_log = Rf_asLogical(log_probability);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        double standardized;

        if (strategy_code == 0) {
            standardized = quantiles[i] - locations[i] + INLACIRC_PI;
            standardized -= 2.0 * INLACIRC_PI *
                            floor(standardized / (2.0 * INLACIRC_PI));
            standardized -= INLACIRC_PI;
        } else {
            standardized = quantiles[i] - locations[i];
        }

        result[i] = plogis(
            INLAcirc_approx_linear(standardized, x_grid, y_grid, grid_length),
            0.0, 1.0, 1, return_log);
    }

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_qvm(SEXP p, SEXP mu, SEXP grid_x, SEXP grid_y)
{
    const int count = Rf_length(p);
    const int grid_length = Rf_length(grid_x);
    const double *probabilities = REAL(p);
    const double *locations = REAL(mu);
    const double *x_grid = REAL(grid_x);
    const double *y_grid = REAL(grid_y);
    const double minimum_logit = y_grid[0];
    const double maximum_logit = y_grid[grid_length - 1];
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        double logit_probability;

        if (ISNAN(probabilities[i])) {
            result[i] = NA_REAL;
            continue;
        }
        if (probabilities[i] <= 0.0) {
            result[i] = locations[i] - INLACIRC_PI;
            continue;
        }
        if (probabilities[i] >= 1.0) {
            result[i] = locations[i] + INLACIRC_PI;
            continue;
        }

        logit_probability = qlogis(probabilities[i], 0.0, 1.0, 1, 0);
        if (logit_probability < minimum_logit) {
            logit_probability = minimum_logit;
        }
        if (logit_probability > maximum_logit) {
            logit_probability = maximum_logit;
        }

        result[i] =
            INLAcirc_approx_linear(logit_probability, y_grid, x_grid,
                                   grid_length) +
            locations[i];
    }

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_rvm(SEXP n, SEXP mu, SEXP grid_x, SEXP grid_y)
{
    const int count = Rf_asInteger(n);
    const int grid_length = Rf_length(grid_x);
    const double *locations = REAL(mu);
    const double *x_grid = REAL(grid_x);
    const double *y_grid = REAL(grid_y);
    const double minimum_logit = y_grid[0];
    const double maximum_logit = y_grid[grid_length - 1];
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    GetRNGstate();
    for (int i = 0; i < count; ++i) {
        double logit_probability = qlogis(unif_rand(), 0.0, 1.0, 1, 0);

        if (logit_probability < minimum_logit) {
            logit_probability = minimum_logit;
        }
        if (logit_probability > maximum_logit) {
            logit_probability = maximum_logit;
        }
        result[i] =
            INLAcirc_approx_linear(logit_probability, y_grid, x_grid,
                                   grid_length) +
            locations[i];
    }
    PutRNGstate();

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_plavm(SEXP q,
                      SEXP eta,
                      SEXP log_probability,
                      SEXP grid_x,
                      SEXP grid_y)
{
    const int count = Rf_length(q);
    const int grid_length = Rf_length(grid_x);
    const double *quantiles = REAL(q);
    const double *predictors = REAL(eta);
    const double *x_grid = REAL(grid_x);
    const double *y_grid = REAL(grid_y);
    const int return_log = Rf_asLogical(log_probability);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        const double shifted = tan(quantiles[i] / 2.0) - predictors[i];
        const double transformed = 2.0 * atan(shifted);
        const double logit_probability =
            INLAcirc_approx_linear(transformed, x_grid, y_grid, grid_length);

        result[i] = plogis(logit_probability, 0.0, 1.0, 1, return_log);
    }

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_qlavm(SEXP p, SEXP eta, SEXP grid_x, SEXP grid_y)
{
    const int count = Rf_length(p);
    const int grid_length = Rf_length(grid_x);
    const double *probabilities = REAL(p);
    const double *predictors = REAL(eta);
    const double *x_grid = REAL(grid_x);
    const double *y_grid = REAL(grid_y);
    const double minimum_logit = y_grid[0];
    const double maximum_logit = y_grid[grid_length - 1];
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    for (int i = 0; i < count; ++i) {
        double probability = probabilities[i];
        double logit_probability;
        double transformed;

        if (ISNAN(probability)) {
            result[i] = NA_REAL;
            continue;
        }
        if (probability <= 0.0) {
            probability = 1e-15;
        }
        if (probability >= 1.0) {
            probability = 1.0 - 1e-15;
        }

        logit_probability = qlogis(probability, 0.0, 1.0, 1, 0);
        if (logit_probability < minimum_logit) {
            logit_probability = minimum_logit;
        }
        if (logit_probability > maximum_logit) {
            logit_probability = maximum_logit;
        }

        transformed = INLAcirc_approx_linear(logit_probability, y_grid,
                                              x_grid, grid_length);
        result[i] =
            2.0 * atan(tan(transformed / 2.0) + predictors[i]);
    }

    UNPROTECT(1);
    return output;
}

SEXP INLAcirc_C_rlavm(SEXP n, SEXP eta, SEXP grid_x, SEXP grid_y)
{
    const int count = Rf_asInteger(n);
    const int grid_length = Rf_length(grid_x);
    const double *predictors = REAL(eta);
    const double *x_grid = REAL(grid_x);
    const double *y_grid = REAL(grid_y);
    const double minimum_logit = y_grid[0];
    const double maximum_logit = y_grid[grid_length - 1];
    SEXP output = PROTECT(Rf_allocVector(REALSXP, count));
    double *result = REAL(output);

    GetRNGstate();
    for (int i = 0; i < count; ++i) {
        double logit_probability = qlogis(unif_rand(), 0.0, 1.0, 1, 0);
        double transformed;

        if (logit_probability < minimum_logit) {
            logit_probability = minimum_logit;
        }
        if (logit_probability > maximum_logit) {
            logit_probability = maximum_logit;
        }

        transformed = INLAcirc_approx_linear(logit_probability, y_grid,
                                              x_grid, grid_length);
        result[i] =
            2.0 * atan(tan(transformed / 2.0) + predictors[i]);
    }
    PutRNGstate();

    UNPROTECT(1);
    return output;
}
