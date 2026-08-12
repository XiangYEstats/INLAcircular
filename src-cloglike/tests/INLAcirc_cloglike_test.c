#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../INLAcirc_cloglike.h"
#include "../INLAcirc_common.h"

static void INLAcirc_test_set_entry(inla_cgeneric_vec_tp *entry,
                                    char *name,
                                    double *value)
{
    entry->name = name;
    entry->len = 1;
    entry->ints = NULL;
    entry->doubles = value;
    entry->chars = NULL;
}

static void INLAcirc_test_close(double actual,
                                double expected,
                                double tolerance)
{
    assert(isfinite(actual));
    assert(fabs(actual - expected) <= tolerance);
}

int main(void)
{
    double link = 0.0;
    double prior = 0.0;
    double prior_u = 0.5;
    double prior_alpha = 0.5;
    double initial = log(5.0);
    double fixed = 0.0;
    double theta[] = {log(5.0)};
    double observation[] = {0.25};
    double predictors[] = {-0.5, 0.0, 0.5};
    double result[3] = {0.0, 0.0, 0.0};
    double *allocated;
    inla_cgeneric_vec_tp entries[6];
    inla_cgeneric_vec_tp *entry_pointers[6];
    inla_cgeneric_data_tp data;

    memset(&data, 0, sizeof(data));
    INLAcirc_test_set_entry(&entries[0], "lavm.link", &link);
    INLAcirc_test_set_entry(&entries[1], "lavm.prior", &prior);
    INLAcirc_test_set_entry(&entries[2], "lavm.u", &prior_u);
    INLAcirc_test_set_entry(&entries[3], "lavm.alpha", &prior_alpha);
    INLAcirc_test_set_entry(&entries[4], "lavm.initial.theta", &initial);
    INLAcirc_test_set_entry(&entries[5], "lavm.fixed.theta", &fixed);
    for (int i = 0; i < 6; ++i) {
        entry_pointers[i] = &entries[i];
    }
    data.n_doubles = 6;
    data.doubles = entry_pointers;

    INLAcirc_test_close(INLAcirc_log_bessel_i0_scaled(5.0),
                        -1.6953182241774667, 1e-12);
    INLAcirc_test_close(
        INLAcirc_log_bessel_i0_scaled_from_log_kappa(log(5.0)),
        -1.6953182241774667, 1e-12);
    INLAcirc_test_close(
        INLAcirc_log_bessel_i0_scaled_from_log_kappa(log(100000.0)),
        -6.6754000156835369, 1e-10);

    allocated = INLAcirc_cloglike_lavm(
        INLA_CLOGLIKE_INITIAL, NULL, &data, 1, observation, 3,
        predictors, result);
    assert(allocated != NULL);
    INLAcirc_test_close(allocated[0], 1.0, 0.0);
    INLAcirc_test_close(allocated[1], log(5.0), 1e-14);
    free(allocated);

    allocated = INLAcirc_cloglike_lavm(
        INLA_CLOGLIKE_LOG_PRIOR, theta, &data, 1, observation, 3,
        predictors, result);
    assert(allocated != NULL);
    assert(isfinite(allocated[0]));
    free(allocated);

    allocated = INLAcirc_cloglike_lavm(
        INLA_CLOGLIKE_LOGLIKE, theta, &data, 1, observation, 3,
        predictors, result);
    assert(allocated == NULL);
    INLAcirc_test_close(result[0], -3.2704584703407269, 1e-8);
    INLAcirc_test_close(result[1], -0.29799673367865465, 1e-8);
    INLAcirc_test_close(result[2], -1.4871405798164146, 1e-8);

    fixed = 1.0;
    allocated = INLAcirc_cloglike_lavm(
        INLA_CLOGLIKE_INITIAL, NULL, &data, 1, observation, 3,
        predictors, result);
    assert(allocated != NULL);
    INLAcirc_test_close(allocated[0], 0.0, 0.0);
    free(allocated);

    allocated = INLAcirc_cloglike_lavm(
        INLA_CLOGLIKE_LOG_PRIOR, NULL, &data, 1, observation, 3,
        predictors, result);
    assert(allocated != NULL);
    INLAcirc_test_close(allocated[0], 0.0, 0.0);
    free(allocated);

    allocated = INLAcirc_cloglike_lavm(
        INLA_CLOGLIKE_LOGLIKE, NULL, &data, 1, observation, 3,
        predictors, result);
    assert(allocated == NULL);
    INLAcirc_test_close(result[0], -3.2704584703407269, 1e-8);
    INLAcirc_test_close(result[1], -0.29799673367865465, 1e-8);
    INLAcirc_test_close(result[2], -1.4871405798164146, 1e-8);

    link = 1.0;
    allocated = INLAcirc_cloglike_lavm(
        INLA_CLOGLIKE_LOGLIKE, NULL, &data, 1, observation, 3,
        predictors, result);
    assert(allocated == NULL);
    for (int i = 0; i < 3; ++i) {
        assert(isfinite(result[i]));
    }

    link = 2.0;
    allocated = INLAcirc_cloglike_lavm(
        INLA_CLOGLIKE_LOGLIKE, NULL, &data, 1, observation, 3,
        predictors, result);
    assert(allocated == NULL);
    for (int i = 0; i < 3; ++i) {
        assert(isfinite(result[i]));
    }

    allocated = INLAcirc_cloglike_lavm(
        INLA_CLOGLIKE_CDF, NULL, &data, 1, observation, 3,
        predictors, result);
    assert(allocated == NULL);
    for (int i = 0; i < 3; ++i) {
        assert(isnan(result[i]));
    }

    puts("INLAcirc cloglike tests passed");
    return 0;
}
