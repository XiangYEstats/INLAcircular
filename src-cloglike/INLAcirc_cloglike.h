#ifndef INLACIRC_CLOGLIKE_H
#       define INLACIRC_CLOGLIKE_H

#if defined(INLACIRC_TEST_CGENERIC_COMPAT)
#include "tests/INLAcirc_cgeneric_test_compat.h"
#else
#include "cgeneric.h"
#endif

#       ifdef __cplusplus
extern "C" {
#       endif

	double *INLAcirc_cloglike_lavm(inla_cloglike_cmd_tp cmd,
				       double *theta, inla_cgeneric_data_tp * data, int ny, double *y, int nx, double *x, double *result);

#       ifdef __cplusplus
}
#       endif
#endif
