suppressPackageStartupMessages(library(INLAcircular))

install_message <- get(
  ".INLAcircular_inla_install_message",
  envir = asNamespace("INLAcircular"),
  inherits = FALSE
)()

stopifnot(
  grepl("not installed", install_message, fixed = TRUE),
  grepl("install.packages(\"INLA\"", install_message, fixed = TRUE),
  grepl("https://inla.r-inla-download.org/R/stable", install_message,
        fixed = TRUE),
  grepl("Attaching INLA is optional", install_message, fixed = TRUE)
)

validate_inla <- get(
  ".INLAcircular_validate_inla",
  envir = asNamespace("INLAcircular"),
  inherits = FALSE
)

capture_error <- function(expression) {
  tryCatch(
    force(expression),
    error = function(error) error
  )
}

old_version_error <- capture_error(validate_inla(
  base::package_version("25.08.20"),
  "inla.cloglike.define"
))
missing_feature_error <- capture_error(validate_inla(
  base::package_version("25.08.21"),
  character()
))

stopifnot(
  isTRUE(validate_inla(
    base::package_version("25.08.21"),
    "inla.cloglike.define"
  )),
  inherits(old_version_error, "error"),
  grepl("requires INLA >= 25.8.21", conditionMessage(old_version_error),
        fixed = TRUE),
  inherits(missing_feature_error, "error"),
  grepl("inla.cloglike.define", conditionMessage(missing_feature_error),
        fixed = TRUE)
)

# Simulate a hard-dependencies-only installation even on a developer machine
# where INLA is installed in a separate user library. INLAcircular is already
# loaded, so temporarily hiding that library tests every public runtime guard
# without uninstalling or modifying INLA.
local({
  original_libraries <- .libPaths()
  on.exit(.libPaths(original_libraries), add = TRUE)

  inla_path <- tryCatch(
    find.package("INLA", quiet = TRUE),
    error = function(e) ""
  )
  if (length(inla_path) == 1L && nzchar(inla_path)) {
    if ("INLA" %in% loadedNamespaces()) {
      unloadNamespace("INLA")
    }
    inla_library <- normalizePath(dirname(inla_path), mustWork = TRUE)
    retained <- original_libraries[
      normalizePath(original_libraries, mustWork = TRUE) != inla_library
    ]
    .libPaths(retained)
  }

  if (!requireNamespace("INLA", quietly = TRUE)) {
    model <- likelihood(
      y ~ intercept(name = "beta0"),
      family = "gaussian"
    )
    errors <- list(
      capture_error(INLAcircular::inla()),
      capture_error(lavm.cloglike()),
      capture_error(inlacc(model, data = data.frame(y = 1))),
      capture_error(summary(structure(list(), class = "inlacc")))
    )

    stopifnot(
      all(vapply(errors, inherits, logical(1L), what = "error")),
      all(vapply(
        errors,
        function(error) identical(conditionMessage(error), install_message),
        logical(1L)
      ))
    )
  }
})
