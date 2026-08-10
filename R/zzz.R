.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion(pkgname)

  # Pull the auto-generated build timestamp
  pkg_date <- utils::packageDescription(pkgname, fields = "Packaged")

  # If the package hasn't been formally built yet, just use today's date
  if (is.null(pkg_date) || is.na(pkg_date)) {
    display_date <- as.character(Sys.Date())
  } else {
    # Extract just the YYYY-MM-DD from the "2026-03-12 11:39:52 UTC; xiang" string
    display_date <- as.character(as.Date(pkg_date))
  }

  inla_version <- tryCatch(
    as.character(utils::packageVersion("INLA")),
    error = function(e) "not installed"
  )

  msg <- paste0(
    "=============================================================\n",
    " INLAcircular: Bayesian Joint Circular Regression with INLA\n",
    " Version: ", pkg_version, " (", display_date, ")\n",
    " INLA version: ", inla_version, "\n",
    "============================================================="
  )

  packageStartupMessage(msg)

  # graphpcor is deliberately optional because it is needed only for the
  # multivariate LKJ covariance model. Check availability without attaching
  # or loading its namespace.
  if (!nzchar(system.file(package = "graphpcor"))) {
    packageStartupMessage(
      "Optional dependency 'graphpcor' is required for model = \"iidkd_LKJ\"."
    )
  }
}
