# INLAcircular

## Installation
You can install the development version of `INLAcircular` from GitHub. 
Because this package depends on `INLA`, ensure your R session knows where to find it:

```R
options(repos = c(getOption("repos"), INLA = "[https://inla.r-inla-download.org/R/testing](https://inla.r-inla-download.org/R/testing)"))
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("XiangYEstats/INLAcircular")
