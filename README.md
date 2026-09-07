# INLAcircular

`INLAcircular` is an R package for Bayesian circular regression and joint circular–linear models using Integrated Nested Laplace Approximation (INLA).

## Installation

After the package is published on CRAN, install its core functionality with:

```r
install.packages("INLAcircular")
```

For a development installation from GitHub, install `remotes` if necessary:

```r
install.packages("remotes")
remotes::install_github(
  "XiangYEstats/INLAcircular"
)
```

The circular distributions and PC-prior functions do not require INLA. Model
fitting uses optional INLA >= 25.08.21 from the current **stable INLA
repository**, which is not installed automatically. Install it separately when
needed:

```r
install.packages(
  "INLA",
  repos = c(
    getOption("repos"),
    INLA = "https://inla.r-inla-download.org/R/stable"
  )
)
```

Attaching INLA is optional because `INLAcircular` uses its namespace directly.
If both packages are attached and you want the compatibility `inla()` function
on the search path, load them in this order:

```r
library(INLA)
library(INLAcircular)
```

## Documentation

[Read the full user guide (PDF)](output/pdf/INLAcircular-guide.pdf)

[Download the HTML guide](output/html/INLAcircular-guide.html?raw=true)
(open the downloaded file in your browser).

The comprehensive guide covers the circular distributions, all PC-prior
densities and d/p/q/r functions, direct fitting with `INLA::inla()`, model
specification, `inlacc()`, and `graphpcor` integration. After installation,
open it with:

```
vignette("INLAcircular-guide", package = "INLAcircular")
```

Function-level documentation is available through
`help(package = "INLAcircular")` and `?function_name`.

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
