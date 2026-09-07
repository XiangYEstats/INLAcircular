#include <R.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rinternals.h>

extern SEXP INLAcirc_C_dvm(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_dlavm(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_dpc_vm0(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_ppc_vm0(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_qpc_vm0(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_rpc_vm0(SEXP, SEXP);
extern SEXP INLAcirc_C_dpc_vminf(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_ppc_vminf(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_qpc_vminf(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_rpc_vminf(SEXP, SEXP);
extern SEXP INLAcirc_C_pvm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_qvm(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_rvm(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_plavm(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_qlavm(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_rlavm(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_bessel_i(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_dcardioid(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_dwrappedcauchy(SEXP, SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_dpc_card0(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_ppc_card0(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_qpc_card0(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_rpc_card0(SEXP, SEXP);
extern SEXP INLAcirc_C_dpc_card(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_ppc_card(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_qpc_card(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_rpc_card(SEXP, SEXP);
extern SEXP INLAcirc_C_dpc_wc(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_ppc_wc(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_qpc_wc(SEXP, SEXP, SEXP);
extern SEXP INLAcirc_C_rpc_wc(SEXP, SEXP);

static const R_CallMethodDef INLAcirc_call_entries[] = {
    {"INLAcirc_C_dvm", (DL_FUNC)&INLAcirc_C_dvm, 4},
    {"INLAcirc_C_dlavm", (DL_FUNC)&INLAcirc_C_dlavm, 4},
    {"INLAcirc_C_dpc_vm0", (DL_FUNC)&INLAcirc_C_dpc_vm0, 3},
    {"INLAcirc_C_ppc_vm0", (DL_FUNC)&INLAcirc_C_ppc_vm0, 3},
    {"INLAcirc_C_qpc_vm0", (DL_FUNC)&INLAcirc_C_qpc_vm0, 3},
    {"INLAcirc_C_rpc_vm0", (DL_FUNC)&INLAcirc_C_rpc_vm0, 2},
    {"INLAcirc_C_dpc_vminf", (DL_FUNC)&INLAcirc_C_dpc_vminf, 3},
    {"INLAcirc_C_ppc_vminf", (DL_FUNC)&INLAcirc_C_ppc_vminf, 3},
    {"INLAcirc_C_qpc_vminf", (DL_FUNC)&INLAcirc_C_qpc_vminf, 3},
    {"INLAcirc_C_rpc_vminf", (DL_FUNC)&INLAcirc_C_rpc_vminf, 2},
    {"INLAcirc_C_pvm", (DL_FUNC)&INLAcirc_C_pvm, 6},
    {"INLAcirc_C_qvm", (DL_FUNC)&INLAcirc_C_qvm, 4},
    {"INLAcirc_C_rvm", (DL_FUNC)&INLAcirc_C_rvm, 4},
    {"INLAcirc_C_plavm", (DL_FUNC)&INLAcirc_C_plavm, 5},
    {"INLAcirc_C_qlavm", (DL_FUNC)&INLAcirc_C_qlavm, 4},
    {"INLAcirc_C_rlavm", (DL_FUNC)&INLAcirc_C_rlavm, 4},
    {"INLAcirc_C_bessel_i", (DL_FUNC)&INLAcirc_C_bessel_i, 3},
    {"INLAcirc_C_dcardioid", (DL_FUNC)&INLAcirc_C_dcardioid, 4},
    {"INLAcirc_C_dwrappedcauchy", (DL_FUNC)&INLAcirc_C_dwrappedcauchy, 4},
    {"INLAcirc_C_dpc_card0", (DL_FUNC)&INLAcirc_C_dpc_card0, 3},
    {"INLAcirc_C_ppc_card0", (DL_FUNC)&INLAcirc_C_ppc_card0, 3},
    {"INLAcirc_C_qpc_card0", (DL_FUNC)&INLAcirc_C_qpc_card0, 3},
    {"INLAcirc_C_rpc_card0", (DL_FUNC)&INLAcirc_C_rpc_card0, 2},
    {"INLAcirc_C_dpc_card", (DL_FUNC)&INLAcirc_C_dpc_card, 3},
    {"INLAcirc_C_ppc_card", (DL_FUNC)&INLAcirc_C_ppc_card, 3},
    {"INLAcirc_C_qpc_card", (DL_FUNC)&INLAcirc_C_qpc_card, 3},
    {"INLAcirc_C_rpc_card", (DL_FUNC)&INLAcirc_C_rpc_card, 2},
    {"INLAcirc_C_dpc_wc", (DL_FUNC)&INLAcirc_C_dpc_wc, 3},
    {"INLAcirc_C_ppc_wc", (DL_FUNC)&INLAcirc_C_ppc_wc, 3},
    {"INLAcirc_C_qpc_wc", (DL_FUNC)&INLAcirc_C_qpc_wc, 3},
    {"INLAcirc_C_rpc_wc", (DL_FUNC)&INLAcirc_C_rpc_wc, 2},
    {NULL, NULL, 0}
};

void attribute_visible R_init_INLAcircular(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, INLAcirc_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
