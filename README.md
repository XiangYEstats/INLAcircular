# INLAcircular <img src="man/figures/logo.png" align="right" width="160" alt="INLAcircular logo" />

`INLAcircular` is an R package for Bayesian circular regression and joint circular-linear models using Integrated Nested Laplace Approximation (INLA). It provides circular distributions, penalized-complexity priors, structured latent effects, shared predictors, and optional LKJ multivariate random effects.

## Documentation

Start with the quick example below, then use the detailed guide for data preparation, priors, structured effects, joint models, LKJ effects, and interpreting results. The guide is the repository's user manual and the source of the installed R vignette.

- [User guide (PDF)](output/pdf/INLAcircular-guide.pdf): installation, complete examples, input and output conventions, priors, troubleshooting, and reproducibility.
- [User guide (HTML)](output/html/INLAcircular-guide.html): on GitHub, download the raw file and open it in a browser.
- [Function reference (PDF)](output/pdf/INLAcircular-reference.pdf): the manual generated from the same help files used by R and CRAN.
- [Guide source](vignettes/INLAcircular-guide.Rmd) and [longer developer examples](tests/manual/README.md).

The repository guide and the installed R vignette share one source. Their model outputs depend on the build mode, which is stated in the guide. The ordinary package build runs the lightweight examples; repository builds can also execute all INLA examples. After installing a package with built vignettes, open the guide with:

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
fit <- inlacc(x_model, data = dat, num.threads = "2:1")
summary(fit)
```

The simulation uses an intercept of `0.5`, a coefficient of `1` for `z`, and concentration `kappa = 50`. The fitted values will vary around these values. The result contains posterior coefficients on the predictor scale and concentration on the natural `kappa` scale.

For a single response, the package also provides an `inla()` interface. Using the same `dat`:

```r
fit_inla <- INLAcircular::inla(
  x ~ z,
  family = "lavm",
  data = dat,
  num.threads = "2:1"
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

The dataset contains 8,759 hourly climate-normal records for JFK International Airport, with wind direction, speed, temperature, and calendar fields. Its timestamps use 2010 as a reference calendar, not an observation year. The guide documents the actual columns and prepares a runnable example, including `hour + 1L` for latent node IDs `1:24`.

For example, construct the hourly mapping for an `f(hour_effect, ...)` term with:

```r
hour_index <- index(
  "hour_effect",
  data.id = wind_newyork$hour + 1L,
  process.id = 1:24
)
```

Pass this mapping as `latent.index = hour_index` when fitting the corresponding model with `inlacc()`. Build `data.id` from the same rows, in the same order, as the fitting data; recreate the mapping after subsetting. The guide includes the complete cyclic model. The original wind-speed and temperature units and the climate-normal baseline period still need to be established from the upstream source before physical interpretation.

## Build the documentation and run checks

GitHub installations skip vignette building by default. To include the installed HTML vignette, install `knitr` and `rmarkdown`, ensure Pandoc is available, and use:

```r
remotes::install_github(
  "XiangYEstats/INLAcircular",
  build_vignettes = TRUE,
  force = TRUE
)
```

To rebuild both repository guides from the current source, also install `pkgload` and a LaTeX distribution with `pdflatex`. From the repository root:

```sh
Rscript tools/render-guide.R
Rscript tools/render-guide.R --run-models
Rscript tools/render-guide.R --run-lkj
```

Choose the first command for the lightweight edition, the second for ordinary INLA examples, or the third for all examples including LKJ. Enabled examples execute separately for HTML and PDF. Missing dependencies for a requested mode cause an error.

The CRAN-style reference manual can be regenerated with:

```sh
R CMD Rd2pdf --no-preview --force --output=output/pdf/INLAcircular-reference.pdf .
```

For standard package checks:

```sh
R CMD build .
R CMD check --as-cran INLAcircular_1.0.0.tar.gz
```

Use the archive filename printed by the build if the version changes. [Automated tests](tests) cover numerical and interface behaviour; [manual scripts](tests/manual/README.md) exercise larger models. [CI configuration](.github/workflows/R-CMD-check.yaml) defines platform checks; consult its run results for the particular version being used.

## Reproducibility, citation, and support

For a publication, use the exact release tag or commit URL for both source code and its guide. Links under `main` can change. Record the software version, any local changes, dependency versions, seeds, and `sessionInfo()`; a versioned archive DOI can be included if one has been created.

```r
citation("INLAcircular")
sessionInfo()
```

The maintained method references are in [inst/CITATION](inst/CITATION). Report reproducible problems through [GitHub issues](https://github.com/XiangYEstats/INLAcircular/issues), including a small example and environment details. Maintainer contact is available with `maintainer("INLAcircular")`.

## License

GNU General Public License, version 3 or later. See [LICENSE.md](LICENSE.md) and the [third-party copyright notices](inst/COPYRIGHTS), including notices for the bundled Bessel routines and previously MIT-licensed original material. The installed package includes the same GPL text under `LICENSES/GPL-3`.
