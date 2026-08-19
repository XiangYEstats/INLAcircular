#ifndef INLACIRC_CGENERIC_COMPAT_H
#define INLACIRC_CGENERIC_COMPAT_H

/*
 * Minimal declarations for the public INLA cloglike ABI.
 *
 * This header belongs only to the R-package shared library.  The separate
 * src-cloglike module includes the cgeneric.h supplied by the INLA source
 * tree when it is compiled into inla-build.
 */

typedef enum {
    INLA_CLOGLIKE_INITIAL = 1,
    INLA_CLOGLIKE_LOG_PRIOR,
    INLA_CLOGLIKE_LOGLIKE,
    INLA_CLOGLIKE_CDF,
    INLA_CLOGLIKE_QUIT
} inla_cloglike_cmd_tp;

typedef struct {
    char *name;
    int nrow;
    int ncol;
    double *x;
} inla_cgeneric_mat_tp;

typedef struct {
    char *name;
    int nrow;
    int ncol;
    int n;
    int *i;
    int *j;
    double *x;
} inla_cgeneric_smat_tp;

typedef struct {
    char *name;
    int len;
    int *ints;
    double *doubles;
    char *chars;
} inla_cgeneric_vec_tp;

typedef struct {
    int max;
    int outer;
    int inner;
} inla_cgeneric_threads_tp;

typedef struct {
    inla_cgeneric_threads_tp threads;
    int n_ints;
    inla_cgeneric_vec_tp **ints;
    int n_doubles;
    inla_cgeneric_vec_tp **doubles;
    int n_chars;
    inla_cgeneric_vec_tp **chars;
    int n_mats;
    inla_cgeneric_mat_tp **mats;
    int n_smats;
    inla_cgeneric_smat_tp **smats;
    int processed;
    void *cache;
    int theta_all_n;
    char **theta_all_names;
} inla_cgeneric_data_tp;

#endif
