#include <R.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rinternals.h>

extern SEXP INLAcirc_C_dvm(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_dlavm(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_pvm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_qvm(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_rvm(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_plavm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_qlavm(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_rlavm(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_bessel_i(SEXP, SEXP, SEXP);

static const R_CallMethodDef INLAcirc_call_entries[] = {
    {"INLAcirc_C_dvm", (DL_FUNC)&INLAcirc_C_dvm, 4},
    {"INLAcirc_C_dlavm", (DL_FUNC)&INLAcirc_C_dlavm, 4},
    {"INLAcirc_C_pvm", (DL_FUNC)&INLAcirc_C_pvm, 6},
    {"INLAcirc_C_qvm", (DL_FUNC)&INLAcirc_C_qvm, 4},
    {"INLAcirc_C_rvm", (DL_FUNC)&INLAcirc_C_rvm, 4},
    {"INLAcirc_C_plavm", (DL_FUNC)&INLAcirc_C_plavm, 5},
    {"INLAcirc_C_qlavm", (DL_FUNC)&INLAcirc_C_qlavm, 4},
    {"INLAcirc_C_rlavm", (DL_FUNC)&INLAcirc_C_rlavm, 4},
    {"INLAcirc_C_bessel_i", (DL_FUNC)&INLAcirc_C_bessel_i, 3},
    {NULL, NULL, 0}
};

void attribute_visible R_init_INLAcircular(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, INLAcirc_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
