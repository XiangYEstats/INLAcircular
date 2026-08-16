# Manual developer tests

These scripts are for package development and regression checking. They are
kept under `tests/manual` so that they are present in the source project but
are not installed as user-facing examples and are not run automatically by
`R CMD check`.

Open any numbered `.R` file and use the RStudio **Source** button, or run it
with `Rscript`. No working-directory setup is required. The wind and
biomechanical scripts locate their data relative to their own file paths.

## Scripts

- `01-direct-inla-lavm.R`: direct `inla(..., family = "lavm")` checks.
- `02-joint-circular-linear.R`: LAvM-Gaussian joint model with an RW2 effect.
- `03-joint-cyclic-circular-linear.R`: cyclic RW2 index-mapping check.
- `04-gamma-regression.R`: ordinary non-LAvM family through `inlacc()`.
- `05-multivariate-joint-model.R`: two LAvM, Gaussian, and Poisson responses.
- `06-wind-new-york.R`: New York wind joint circular-Gamma application.
- `07-biomechanical-lkj.R`: six-response biomechanical `iidkd_LKJ` model.

The first and fourth scripts are the quickest smoke tests. The joint, wind,
and biomechanical fits can take substantially longer.

## Data

- `data/wind_newyork_2010.csv` is the cleaned one-year JFK dataset with only
  the eight variables required by the test. Wind direction is in radians.
- `data/biomechanical/` contains the six CSV files required by the LKJ test.

The original session-history file, duplicate/incomplete wind CSV, and source
PDF documents are intentionally omitted because they are not needed to test
the package.
