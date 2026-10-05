# INLAcircular 1.0.0

## Initial CRAN release

- Shows the installed version, automatic build date, and user-guide links
  when the package is attached.
- Preserves missing values in circular CDFs without reading outside the
  interpolation grid.
- Supports a vector sample-size argument in rvm() and rlavm() by using its
  length, as documented.
- Limits regular INLA fitting tests to two threads and keeps developer-only
  manual scripts and duplicate datasets out of the CRAN source archive.
- Records method citations and the attribution and reuse notices for the
  Acklam approximation and INLA interface declarations.
- Provides von Mises and link-adjusted von Mises d/p/q/r functions.
- Provides compiled cardioid and wrapped Cauchy densities.
- Provides a modified Bessel I implementation adapted from GSL 2.8, without
  requiring a system GSL installation. The package is distributed under GPL
  version 3 or later; upstream notices and the earlier MIT permission notice
  for original material are retained.
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
