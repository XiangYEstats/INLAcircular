#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "INLAcirc_cloglike.h"
#include "INLAcirc_common.h"

static double *INLAcirc_alloc_doubles(size_t count)
{
    return (double *)malloc(count * sizeof(double));
}

/* Analytical PC-prior rate calculations; no global cache, hence thread-safe. */
static double INLAcirc_lavm_pc_inf_lambda(double u, double alpha)
{
    return -log(1.0 - alpha) / sqrt(1.0 - u);
}

static double INLAcirc_lavm_pc_0_lambda(double u, double alpha)
{
    double kappa;

    if (u <= 1e-5) {
        return -log(alpha) / sqrt(u);
    }
    if (u < 0.53) {
        kappa = 2.0 * u + u * u * u + 5.0 * pow(u, 5.0) / 6.0;
    } else if (u < 0.85) {
        kappa = -0.4 + 1.39 * u + 0.43 / (1.0 - u);
    } else {
        kappa = 1.0 / (2.0 * (1.0 - u));
    }

    for (int iteration = 0; iteration < 7; ++iteration) {
        const double i0 = INLAcirc_bessel_i(kappa, 0.0, 1);
        const double i1 = INLAcirc_bessel_i(kappa, 1.0, 1);
        const double ratio = i1 / i0;
        const double difference = ratio - u;
        const double derivative = 1.0 - ratio * ratio - ratio / kappa;

        if (fabs(difference) < 1e-8) {
            break;
        }
        kappa -= difference / derivative;
    }

    {
        const double i0_scaled = INLAcirc_bessel_i(kappa, 0.0, 1);
        const double distance_squared =
            kappa * u - kappa - log(i0_scaled);
        const double distance = sqrt(fmax(1e-16, distance_squared));
        return -log(alpha) / distance;
    }
}

static double INLAcirc_lavm_inverse_tangent(double observation,
                                             double eta,
                                             double log_kappa)
{
    double x = observation;
    double tangent_half_x;
    double shifted;
    double kernel = 0.0;

    if (x > INLACIRC_PI - 1e-12) {
        x = INLACIRC_PI - 1e-12;
    }
    if (x < -INLACIRC_PI + 1e-12) {
        x = -INLACIRC_PI + 1e-12;
    }

    tangent_half_x = tan(x / 2.0);
    shifted = tangent_half_x - eta;
    if (shifted * shifted > 0.0) {
        kernel = -2.0 * exp(log_kappa + log(shifted * shifted) -
                            log1p(shifted * shifted));
    }

    return kernel - INLACIRC_LOG_2PI -
           INLAcirc_log_bessel_i0_scaled_from_log_kappa(log_kappa) +
           log1p(tangent_half_x * tangent_half_x) -
           log1p(shifted * shifted);
}

static double INLAcirc_lavm_scaled_logit(double observation,
                                         double eta,
                                         double log_kappa)
{
    const double normalized_x =
        (observation + INLACIRC_PI) / (2.0 * INLACIRC_PI);
    const double shifted = INLAcirc_qlogis(normalized_x) - eta;
    const double probability = INLAcirc_plogis(shifted);
    const double sine_squared =
        pow(sin((2.0 * INLACIRC_PI * probability - INLACIRC_PI) / 2.0),
            2.0);

    return -2.0 * exp(log_kappa + log(sine_squared)) -
           INLACIRC_LOG_2PI -
           INLAcirc_log_bessel_i0_scaled_from_log_kappa(log_kappa) +
           log(probability * (1.0 - probability)) -
           log(normalized_x * (1.0 - normalized_x));
}

static double INLAcirc_lavm_scaled_probit(double observation,
                                          double eta,
                                          double log_kappa)
{
    const double normalized_x =
        (observation + INLACIRC_PI) / (2.0 * INLACIRC_PI);
    const double transformed_x = INLAcirc_qnorm(normalized_x);
    const double shifted = transformed_x - eta;
    const double probability = INLAcirc_pnorm(shifted);
    const double sine_squared =
        pow(sin((2.0 * INLACIRC_PI * probability - INLACIRC_PI) / 2.0),
            2.0);

    return -2.0 * exp(log_kappa + log(sine_squared)) -
           INLACIRC_LOG_2PI -
           INLAcirc_log_bessel_i0_scaled_from_log_kappa(log_kappa) +
           0.5 * (transformed_x * transformed_x - shifted * shifted);
}

static double INLAcirc_lavm_pc_prior_inf(double log_kappa, double lambda)
{
    double distance;
    double log_jacobian;

    if (log_kappa <= 11.512925) {
        const double kappa = exp(log_kappa);
        const double i0 = INLAcirc_bessel_i(kappa, 0.0, 1);
        const double i1 = INLAcirc_bessel_i(kappa, 1.0, 1);
        const double i2 = INLAcirc_bessel_i(kappa, 2.0, 1);
        const double ratio = i1 / i0;

        distance = sqrt(fmax(1e-16, 1.0 - ratio));
        log_jacobian = -INLACIRC_LOG_2 - log(distance) +
                       log(fmax(1e-16,
                                (i0 + i2) / (2.0 * i0) - ratio * ratio));
    } else {
        const double inverse_kappa = exp(-log_kappa);
        const double inverse_kappa_squared = inverse_kappa * inverse_kappa;
        const double log_distance =
            -INLACIRC_LOG_2 / 2.0 - log_kappa / 2.0 +
            0.125 * inverse_kappa +
            (7.0 / 64.0) * inverse_kappa_squared;

        distance = exp(log_distance);
        log_jacobian =
            -INLACIRC_LOG_2 - log_distance +
            (-INLACIRC_LOG_2 - 2.0 * log_kappa +
             0.5 * inverse_kappa + 0.625 * inverse_kappa_squared);
    }

    return log(lambda) - lambda * distance + log_jacobian + log_kappa;
}

static double INLAcirc_lavm_pc_prior_0(double log_kappa, double lambda)
{
    const double kappa = exp(log_kappa);
    double distance;
    double log_jacobian_partial;

    if (kappa < 1e-4) {
        double density =
            lambda / 2.0 - (lambda * lambda * kappa) / 4.0 +
            (-9.0 * lambda / 64.0 + lambda * lambda * lambda / 16.0) *
                kappa * kappa;
        if (density <= 1e-16) {
            density = 1e-16;
        }
        return log(density) + log_kappa;
    }

    if (log_kappa <= 11.512925) {
        const double i0 = INLAcirc_bessel_i(kappa, 0.0, 1);
        const double i1 = INLAcirc_bessel_i(kappa, 1.0, 1);
        const double i2 = INLAcirc_bessel_i(kappa, 2.0, 1);
        const double ratio = i1 / i0;
        const double distance_squared =
            kappa * ratio - kappa - log(i0);

        distance = sqrt(fmax(1e-16, distance_squared));
        log_jacobian_partial =
            log(kappa) +
            log(fmax(1e-16,
                     (i0 + i2) / (2.0 * i0) - ratio * ratio));
    } else {
        const double inverse_kappa = exp(-log_kappa);
        const double inverse_kappa_squared = inverse_kappa * inverse_kappa;
        const double distance_squared =
            0.5 * INLACIRC_LOG_2PI - 0.5 + 0.5 * log_kappa -
            0.25 * inverse_kappa - 0.1875 * inverse_kappa_squared -
            (25.0 / 96.0) * inverse_kappa_squared * inverse_kappa -
            (65.0 / 128.0) * inverse_kappa_squared * inverse_kappa_squared -
            (3219.0 / 2560.0) * inverse_kappa_squared *
                inverse_kappa_squared * inverse_kappa;

        distance = sqrt(fmax(1e-16, distance_squared));
        log_jacobian_partial =
            -INLACIRC_LOG_2 - log_kappa + 0.5 * inverse_kappa +
            0.625 * inverse_kappa_squared +
            (59.0 / 48.0) * inverse_kappa_squared * inverse_kappa +
            (203.0 / 64.0) * inverse_kappa_squared *
                inverse_kappa_squared +
            (12743.0 / 1280.0) * inverse_kappa_squared *
                inverse_kappa_squared * inverse_kappa;
    }

    return log(lambda) - lambda * distance - INLACIRC_LOG_2 -
           log(distance) + log_jacobian_partial + log_kappa;
}

double *INLAcirc_cloglike_lavm(inla_cloglike_cmd_tp cmd,
                               double *theta,
                               inla_cgeneric_data_tp *data,
                               int ny,
                               double *y,
                               int nx,
                               double *x,
                               double *result)
{
    double *return_value = NULL;
    double initial_log_kappa = 6.0;
    double prior_u = 0.5;
    double prior_alpha = 0.5;
    double log_kappa;
    int link_code = 0;
    int prior_code = 0;
    int fixed_log_kappa = 0;

    (void)ny;

    if (data != NULL) {
        for (int i = 0; i < data->n_doubles; ++i) {
            const inla_cgeneric_vec_tp *entry = data->doubles[i];

            if (entry == NULL || entry->len < 1 || entry->doubles == NULL) {
                continue;
            }
            if (strcmp(entry->name, "lavm.link") == 0) {
                link_code = (int)entry->doubles[0];
            } else if (strcmp(entry->name, "lavm.prior") == 0) {
                prior_code = (int)entry->doubles[0];
            } else if (strcmp(entry->name, "lavm.u") == 0) {
                prior_u = entry->doubles[0];
            } else if (strcmp(entry->name, "lavm.alpha") == 0) {
                prior_alpha = entry->doubles[0];
            } else if (strcmp(entry->name, "lavm.initial.theta") == 0) {
                initial_log_kappa = entry->doubles[0];
            } else if (strcmp(entry->name, "lavm.fixed.theta") == 0) {
                fixed_log_kappa = (int)entry->doubles[0];
            }
        }
    }

    log_kappa = fixed_log_kappa ? initial_log_kappa : 0.0;
    if (!fixed_log_kappa && theta != NULL) {
        log_kappa = theta[0];
    }
    if (log_kappa > 30.0) {
        log_kappa = 30.0;
    }
    if (log_kappa < -30.0) {
        log_kappa = -30.0;
    }

    switch (cmd) {
    case INLA_CLOGLIKE_INITIAL:
        if (fixed_log_kappa) {
            return_value = INLAcirc_alloc_doubles(1);
            if (return_value != NULL) {
                return_value[0] = 0.0;
            }
        } else {
            return_value = INLAcirc_alloc_doubles(2);
            if (return_value != NULL) {
                return_value[0] = 1.0;
                return_value[1] = initial_log_kappa;
            }
        }
        break;

    case INLA_CLOGLIKE_LOG_PRIOR:
        return_value = INLAcirc_alloc_doubles(1);
        if (return_value != NULL) {
            if (fixed_log_kappa || theta == NULL) {
                return_value[0] = 0.0;
            } else {
                const double lambda =
                    (prior_code == 1)
                        ? INLAcirc_lavm_pc_0_lambda(prior_u, prior_alpha)
                        : INLAcirc_lavm_pc_inf_lambda(prior_u, prior_alpha);
                return_value[0] =
                    (prior_code == 1)
                        ? INLAcirc_lavm_pc_prior_0(log_kappa, lambda)
                        : INLAcirc_lavm_pc_prior_inf(log_kappa, lambda);
            }
        }
        break;

    case INLA_CLOGLIKE_LOGLIKE:
        if ((!fixed_log_kappa && theta == NULL) || y[0] != y[0]) {
            for (int i = 0; i < nx; ++i) {
                result[i] = 0.0;
            }
            break;
        }
        for (int i = 0; i < nx; ++i) {
            if (link_code == 0) {
                result[i] = INLAcirc_lavm_inverse_tangent(
                    y[0], x[i], log_kappa);
            } else if (link_code == 1) {
                result[i] = INLAcirc_lavm_scaled_logit(
                    y[0], x[i], log_kappa);
            } else {
                result[i] = INLAcirc_lavm_scaled_probit(
                    y[0], x[i], log_kappa);
            }
        }
        break;

    case INLA_CLOGLIKE_CDF:
        for (int i = 0; i < nx; ++i) {
            result[i] = NAN;
        }
        break;

    case INLA_CLOGLIKE_QUIT:
        break;
    }

    return return_value;
}
