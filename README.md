# INLAcircular

`INLAcircular` is an R package for Bayesian circular regression and joint circular–linear models using Integrated Nested Laplace Approximation (INLA).

## Installation

`INLAcircular` uses the **INLA testing version**. Install it first:

```
install.packages(
  "INLA",
  repos = c(
    getOption("repos"),
    INLA = "https://inla.r-inla-download.org/R/testing"
  )
)
```

Install `remotes` if necessary:

```
install.packages("remotes")
```

Then install `INLAcircular` from GitHub:

```
remotes::install_github(
  "XiangYEstats/INLAcircular",
  dependencies = TRUE
)
```

Load the package:

```
library(INLAcircular)
```

For models using the LKJ multivariate random-effect model, install `graphpcor`:

```
install.packages("graphpcor")
```

## New York wind dataset

The package includes the processed full-year New York wind dataset:

```
data("wind_newyork", package = "INLAcircular")

head(wind_newyork)
```

The dataset contains wind direction in radians, wind speed, temperature, and time variables for the complete year.

## LKJ multivariate random effects

`INLAcircular` supports an LKJ-based multivariate IID random-effect model through `graphpcor`:

```
f(i, model = "iidkd_LKJ")
```

Its default prior settings are:

```
pc.prior.u     = 1
pc.prior.alpha = 0.5
LKJ.eta        = 5
```

You may override them when needed:

```
f(
  i,
  model = "iidkd_LKJ",
  pc.prior.u = 0.5,
  pc.prior.alpha = 0.01
)
```

and supply `LKJ.eta` in `inlacc()`:

```
fit <- inlacc(
  model,
  data = data,
  LKJ.eta = 2
)
```

## License

MIT License.
