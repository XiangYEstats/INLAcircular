# INLAcircular

`INLAcircular` is an R package for Bayesian circular regression and joint circular–linear models using Integrated Nested Laplace Approximation (INLA).

## Installation

`INLAcircular` uses the **INLA testing version**. Install it first:

```r
install.packages(
  "INLA",
  repos = c(
    getOption("repos"),
    INLA = "https://inla.r-inla-download.org/R/testing"
  )
)
```

Install `remotes` if necessary:

```r
install.packages("remotes")
```

Then install `INLAcircular` from GitHub:

```r
remotes::install_github(
  "XiangYEstats/INLAcircular",
  dependencies = TRUE
)
```

Load the package:

```r
library(INLAcircular)
```

For models using the LKJ multivariate random-effect model, install `graphpcor`:

```r
install.packages("graphpcor")
```

## Included examples

List all installed example scripts:

```r
list.files(
  system.file("examples", package = "INLAcircular"),
  pattern = "\\.R$"
)
```

View an example without running it:

```r
file.show(system.file(
  "examples",
  "wind_new_york_joint_model.R",
  package = "INLAcircular"
))
```

Run the New York wind-data example:

```r
source(system.file(
  "examples",
  "wind_new_york_joint_model.R",
  package = "INLAcircular",
  mustWork = TRUE
))
```

Run a simulated joint circular–linear example:

```r
source(system.file(
  "examples",
  "simulated_joint_circular_linear.R",
  package = "INLAcircular",
  mustWork = TRUE
))
```

## New York wind dataset

The package includes the processed full-year New York wind dataset:

```r
data("wind_newyork", package = "INLAcircular")

head(wind_newyork)
```

The dataset contains wind direction in radians, wind speed, temperature, and time variables for the complete year.

## LKJ multivariate random effects

`INLAcircular` supports an LKJ-based multivariate IID random-effect model through `graphpcor`:

```r
f(i, model = "iidkd_LKJ")
```

Its default prior settings are:

```r
pc.prior.u     = 1
pc.prior.alpha = 0.5
LKJ.eta        = 5
```

You may override them when needed:

```r
f(
  i,
  model = "iidkd_LKJ",
  pc.prior.u = 0.5,
  pc.prior.alpha = 0.01
)
```

and supply `LKJ.eta` in `inlacc()`:

```r
fit <- inlacc(
  model,
  data = data,
  LKJ.eta = 2
)
```

## License

MIT License.
