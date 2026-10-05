#ifndef INLACIRC_CGENERIC_COMPAT_H
#define INLACIRC_CGENERIC_COMPAT_H

/*
 * Minimal declarations for the public INLA cloglike ABI.
 * Adapted from inlaprog/src/cgeneric.h in hrue/r-inla, verified at commit
 * 6afa911f01e975c4e6c921f22288e80d8e70c88f. See inst/COPYRIGHTS for provenance.
 *
 * Copyright 2025 Haavard Rue
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
