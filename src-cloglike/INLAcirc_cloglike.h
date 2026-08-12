#ifndef INLACIRC_CLOGLIKE_H
#define INLACIRC_CLOGLIKE_H

#include "INLAcirc_cgeneric_compat.h"

#ifdef __cplusplus
extern "C" {
#endif

double *INLAcirc_cloglike_lavm(inla_cloglike_cmd_tp cmd,
                               double *theta,
                               inla_cgeneric_data_tp *data,
                               int ny,
                               double *y,
                               int nx,
                               double *x,
                               double *result);

#ifdef __cplusplus
}
#endif

#endif
