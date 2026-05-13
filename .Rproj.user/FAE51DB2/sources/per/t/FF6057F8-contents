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

  inla_version <- utils::packageVersion("INLA")

  msg <- paste0(
    "=============================================================\n",
    " INLAcircular: Bayesian Joint Circular Regression with INLA\n",
    " Version: ", pkg_version, " (", display_date, ")\n",
    " Currently running on INLA version: ", inla_version, "\n",
    "============================================================="
  )

  packageStartupMessage(msg)
}
