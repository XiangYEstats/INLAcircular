#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern SEXP C_dvm(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_dlavm(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_pvm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_qvm(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_rvm(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_plavm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_qlavm(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_rlavm(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_bessel_i(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"C_dvm",      (DL_FUNC) &C_dvm,      4},
  {"C_dlavm",    (DL_FUNC) &C_dlavm,    4},
  {"C_pvm",      (DL_FUNC) &C_pvm,      6},
  {"C_qvm",      (DL_FUNC) &C_qvm,      4},
  {"C_rvm",      (DL_FUNC) &C_rvm,      4},
  {"C_plavm",    (DL_FUNC) &C_plavm,    5},
  {"C_qlavm",    (DL_FUNC) &C_qlavm,    4},
  {"C_rlavm",    (DL_FUNC) &C_rlavm,    4},
  {"C_bessel_i", (DL_FUNC) &C_bessel_i, 3},
  {NULL, NULL, 0}
};

void R_init_INLAcircular(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
