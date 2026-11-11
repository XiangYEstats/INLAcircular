# INLAcircular <img src="man/figures/logo.png" align="right" width="160" alt="INLAcircular logo" />

`INLAcircular` is an R package for Bayesian circular regression and joint circular-linear models using Integrated Nested Laplace Approximation (INLA). It provides circular distributions, penalized-complexity priors, structured latent effects, shared predictors, and optional LKJ multivariate random effects.

## Documentation

Start with the quick example below, then use the detailed guide for data preparation, priors, structured effects, joint models, LKJ effects, and interpreting results.

- [User guide (PDF)](output/pdf/INLAcircular-guide.pdf): installation, complete examples, input and output conventions, priors, troubleshooting, and reproducibility.
- [User guide (HTML)](output/html/INLAcircular-guide.html): on GitHub, download the raw file and open it in a browser.
- [Function reference (PDF)](output/pdf/INLAcircular-reference.pdf): arguments, return values, and examples for individual functions.

After installing from CRAN, open the vignette and function help in R:

```r
vignette("INLAcircular-guide", package = "INLAcircular")
help(package = "INLAcircular")
```

## Installation

The core package requires R >= 4.1.0. Install the released version from [CRAN](https://cran.r-project.org/web/packages/INLAcircular/index.html):

```r
install.packages("INLAcircular")
```

For the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("XiangYEstats/INLAcircular")
```

The numerical functions do not require INLA. Model fitting requires optional INLA >= 25.08.21, installed separately from its stable repository:

```r
install.packages(
  "INLA",
  repos = c(
    getOption("repos"),
    INLA = "https://inla.r-inla-download.org/R/stable"
  )
)
library(INLAcircular)
```

Attaching INLA is optional. If attaching both packages, load `INLA` first and `INLAcircular` second; alternatively use `INLAcircular::inla()` explicitly. The public circular fitting interface is `family = "lavm"`.

LKJ multivariate effects additionally require:

```r
install.packages(c("graphpcor", "INLAtools"))
```

Source installations compile C code and require R's platform-specific build tools. The core package does not require an external GSL installation. Optional packages have their own platform and R-version requirements; see the guide for details.

## Quick start

Check the numerical API without INLA:

```r
library(INLAcircular)
dvm(0, mu = 0, kappa = 0)  # 1 / (2 * pi)
ppc.vm0(qpc.vm0(0.5, lambda = 2), lambda = 2)  # approximately 0.5
```

With INLA installed, run a complete circular regression:

```r
library(INLAcircular)
set.seed(1)
z <- rnorm(60)
x <- rlavm(60, eta = 0.5 + z, kappa = 50)
dat <- data.frame(x = x, z = z)
x_model <- likelihood(x ~ z, family = "lavm")
fit <- inlacc(x_model, data = dat)
summary(fit)
```

The simulation uses an intercept of `0.5`, a coefficient of `1` for `z`, and concentration `kappa = 50`. The fitted values will vary around these values. The result contains posterior coefficients on the predictor scale and concentration on the natural `kappa` scale.

For a single response, the package also provides an `inla()` interface. Using the same `dat`:

```r
fit_inla <- INLAcircular::inla(
  x ~ z,
  family = "lavm",
  data = dat
)
fit_inla$summary.fixed
fit_inla$summary.hyperpar
```

`INLAcircular::inla()` selects the package's compatibility layer explicitly and prepares the circular response automatically. Use `inlacc()` to assemble likelihood blocks, copied predictors, and joint models as shown in the guide.

Supply circular observations as `x` in radians inside `(-pi, pi)`. With the default inverse-tangent link, `2 * atan(eta)` is the circular location parameter corresponding to linear predictor `eta`. See the guide for wrapping, reference directions, and uncertainty summaries.

## Included wind data

```r
data("wind_newyork", package = "INLAcircular")
dim(wind_newyork)
```

The dataset contains 8,759 hourly climate-normal records for JFK International Airport, with wind direction, speed, temperature, and calendar fields. Its timestamps use 2010 as a reference calendar. The guide explains the data's limitations and provides a complete circular regression example with cyclic hourly effects, using `hour + 1L` and `process.id = 1:24`.

## Reproducibility, citation, and support

Use the package citation and record your software environment when reporting an analysis:

```r
citation("INLAcircular")
sessionInfo()
```

For development versions, also record the Git commit used. Report problems through [GitHub issues](https://github.com/XiangYEstats/INLAcircular/issues), including a small reproducible example and the output of `sessionInfo()`. Maintainer contact is available with `maintainer("INLAcircular")`.

## License

GNU General Public License, version 3 or later. See [LICENSE.md](LICENSE.md) and the [third-party copyright notices](inst/COPYRIGHTS), including notices for the bundled Bessel routines and previously MIT-licensed original material. The installed package includes the same GPL text under `LICENSES/GPL-3`.
