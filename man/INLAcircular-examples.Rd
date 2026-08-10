\name{INLAcircular-examples}
\alias{INLAcircular-examples}
\title{Installed INLAcircular Examples}
\description{
Complete scripts illustrating simulated and applied models fitted with
\pkg{INLAcircular}.
}
\details{
The package installs the following scripts under \code{examples/}:
\itemize{
  \item \code{simulated_joint_circular_linear.R}: joint LAvM-Gaussian
  regression with an RW2 effect.
  \item \code{simulated_cyclic_joint_circular_linear.R}: joint regression
  with a cyclic RW2 effect and explicit index mapping.
  \item \code{simulated_gamma_regression.R}: single-response Gamma
  regression.
  \item \code{simulated_multivariate_joint_model.R}: a joint model for two
  circular, one Gaussian, and one Poisson response.
  \item \code{wind_new_york_joint_model.R}: joint circular-Gamma regression
  for New York hourly wind normals, using the processed full-year
  \code{wind_newyork} package dataset.
}
}
\examples{
\dontrun{
source(system.file(
  "examples",
  "simulated_joint_circular_linear.R",
  package = "INLAcircular",
  mustWork = TRUE
))
}
}
\keyword{models}
