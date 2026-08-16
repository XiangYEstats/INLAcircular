# INLAcircular examples

The installed example scripts are:

- `direct_inla_lavm.R`: direct single-response LAvM regression through
  `inla()`, with both default and explicit `control.family` settings.
- `simulated_joint_circular_linear.R`: a joint LAvM-Gaussian model with a
  shared latent predictor and an RW2 effect.
- `simulated_cyclic_joint_circular_linear.R`: the same joint structure with a
  cyclic RW2 effect and explicit index mapping.
- `simulated_gamma_regression.R`: a single-response Gamma regression.
- `simulated_multivariate_joint_model.R`: two circular responses, one Gaussian
  response, and one Poisson response in a joint model.
- `wind_new_york_joint_model.R`: a joint circular-Gamma application to New York
  hourly wind normals. The script loads the processed full-year
  `wind_newyork` package dataset and selects the first quarter for fitting.

For example, run the direct `inla()` script after installing the package with:

```r
source(system.file(
  "examples",
  "direct_inla_lavm.R",
  package = "INLAcircular",
  mustWork = TRUE
))
```
