.INLAcircular_inla_repository <-
  "https://inla.r-inla-download.org/R/stable"

.INLAcircular_minimum_inla_version <- base::package_version("25.08.21")

.INLAcircular_inla_install_message <- function() {
  paste0(
    "Package 'INLA' is required for model fitting but is not installed.\n",
    "Install the stable INLA release first with:\n",
    "  install.packages(\"INLA\", repos = c(getOption(\"repos\"),\n",
    "    INLA = \"", .INLAcircular_inla_repository, "\"))\n",
    "Attaching INLA is optional; INLAcircular loads its namespace when needed."
  )
}

.INLAcircular_inla_upgrade_message <- function(installed_version) {
  paste0(
    "The installed 'INLA' version (", installed_version,
    ") is not compatible with INLAcircular.\n",
    "INLAcircular requires INLA >= ",
    as.character(.INLAcircular_minimum_inla_version),
    " with the exported function 'inla.cloglike.define'.\n",
    "Install or update the stable INLA release with:\n",
    "  install.packages(\"INLA\", repos = c(getOption(\"repos\"),\n",
    "    INLA = \"", .INLAcircular_inla_repository, "\"))"
  )
}

.INLAcircular_validate_inla <- function(version, exports) {
  compatible_version <- tryCatch(
    base::package_version(as.character(version)) >=
      .INLAcircular_minimum_inla_version,
    error = function(error) FALSE
  )
  has_cloglike <- is.character(exports) &&
    "inla.cloglike.define" %in% exports

  if (!isTRUE(compatible_version) || !has_cloglike) {
    stop(
      .INLAcircular_inla_upgrade_message(as.character(version)),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.INLAcircular_require_inla <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(.INLAcircular_inla_install_message(), call. = FALSE)
  }

  .INLAcircular_validate_inla(
    utils::packageVersion("INLA"),
    getNamespaceExports("INLA")
  )

  invisible(TRUE)
}

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

  if (identical(inla_version, "not installed")) {
    packageStartupMessage(.INLAcircular_inla_install_message())
  } else if (base::package_version(inla_version) <
             .INLAcircular_minimum_inla_version) {
    packageStartupMessage(.INLAcircular_inla_upgrade_message(inla_version))
  }

  # graphpcor is deliberately optional because it is needed only for the
  # multivariate LKJ covariance model. Check availability without attaching
  # or loading its namespace.
  if (!nzchar(system.file(package = "graphpcor"))) {
    packageStartupMessage(
      "Optional dependency 'graphpcor' is required for model = \"iidkd_LKJ\"."
    )
  }
}
