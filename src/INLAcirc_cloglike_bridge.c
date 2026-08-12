/*
 * Backward-compatibility build bridge.
 *
 * R builds only translation units found in src/.  Including the independent
 * module here keeps lavm.cloglike() usable from the installed package while
 * preserving src-cloglike/ as the canonical, R-free source for inla-build.
 * INLA must compile src-cloglike/INLAcirc_cloglike.c directly and must not
 * compile this bridge.
 */
#include "../src-cloglike/INLAcirc_cloglike.c"
