# INLAcircular 1.0.0

## Initial CRAN release

- Provides von Mises and link-adjusted von Mises d/p/q/r functions.
- Provides compiled cardioid and wrapped Cauchy densities.
- Provides d/p/q/r PC priors for von Mises, cardioid, and wrapped Cauchy
  concentration parameters.
- Provides direct R-INLA and `inlacc()` interfaces for LAvM and joint models.
- Makes response-predictor copying automatic in `covariate(response)` while
  retaining explicit `predictor = TRUE` and observed-value
  `predictor = FALSE` modes.
- Provides optional graphpcor-backed LKJ multivariate random effects.
- Adds an installed HTML vignette plus automatically knitted repository HTML
  and PDF versions of the comprehensive user guide and function reference.
- Keeps INLA optional at installation time, requires INLA >= 25.08.21 only
  for model fitting, and gives actionable install/upgrade instructions when
  that support is unavailable or incompatible.
- Tests the native mathematical constants and likelihood code with both GCC
  and Clang as part of the cross-platform release workflow.
