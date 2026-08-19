# End-to-end test for the standalone LAvM cloglike source.
#
# Run from any directory with:
#   Rscript src-cloglike/tests/test-cloglike-inla.R
#
# The script compiles src-cloglike/INLAcirc_cloglike.c into a temporary
# shared library, supplies that library to inla.cloglike.define(), and fits a
# small model with INLA.  It does not use the package DLL.

if (!requireNamespace("INLA", quietly = TRUE)) {
  stop("This test requires the 'INLA' package.", call. = FALSE)
}

.INLAcirc_find_test_file <- function() {
  filename <- "test-cloglike-inla.R"
  candidates <- character()

  for (frame in rev(sys.frames())) {
    path <- frame$ofile
    if (is.character(path) && length(path) == 1L && nzchar(path)) {
      candidates <- c(candidates, path)
    }
  }

  file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_argument)) {
    candidates <- c(
      candidates,
      sub("^--file=", "", file_argument[[1L]])
    )
  }

  source_root <- Sys.getenv("INLACIRC_SOURCE_ROOT", unset = "")
  candidates <- c(
    candidates,
    file.path(getwd(), filename),
    file.path(getwd(), "src-cloglike", "tests", filename),
    if (nzchar(source_root)) {
      file.path(source_root, "src-cloglike", "tests", filename)
    }
  )

  candidates <- candidates[
    !is.na(candidates) & nzchar(candidates) & file.exists(candidates)
  ]
  if (!length(candidates)) {
    stop(
      paste0(
        "Cannot locate ", filename, ". Run the complete file with Rscript ",
        "or set INLACIRC_SOURCE_ROOT to the package source directory."
      ),
      call. = FALSE
    )
  }

  normalizePath(candidates[[1L]], mustWork = TRUE)
}

.INLAcirc_run_test <- function() {
  test_file <- .INLAcirc_find_test_file()
  test_directory <- dirname(test_file)
  source_directory <- normalizePath(
    file.path(test_directory, ".."),
    mustWork = TRUE
  )

  build_directory <- tempfile("INLAcirc-cloglike-test-")
  dir.create(build_directory)
  on.exit(unlink(build_directory, recursive = TRUE, force = TRUE), add = TRUE)

  source_files <- c(
    "INLAcirc_cloglike.c",
    "INLAcirc_cloglike.h",
    "INLAcirc_common.h"
  )
  dir.create(file.path(build_directory, "tests"))
  copied <- c(
    file.copy(
      file.path(source_directory, source_files),
      file.path(build_directory, source_files),
      overwrite = TRUE
    ),
    file.copy(
      file.path(test_directory, "INLAcirc_cgeneric_test_compat.h"),
      file.path(
        build_directory,
        "tests",
        "INLAcirc_cgeneric_test_compat.h"
      ),
      overwrite = TRUE
    )
  )
  if (!all(copied)) {
    stop("Could not stage the standalone cloglike sources.", call. = FALSE)
  }

  shared_library <- file.path(
    build_directory,
    paste0("INLAcirc_cloglike", .Platform$dynlib.ext)
  )
  source_file <- file.path(build_directory, "INLAcirc_cloglike.c")
  r_executable <- file.path(R.home("bin"), "R")

  old_working_directory <- setwd(build_directory)
  compile_output <- tryCatch(
    system2(
      r_executable,
      c(
        "CMD", "SHLIB", shQuote(source_file),
        "-o", shQuote(shared_library)
      ),
      stdout = TRUE,
      stderr = TRUE,
      env = c(
        "PKG_CPPFLAGS=-DINLACIRC_TEST_CGENERIC_COMPAT",
        "PKG_LIBS=-lm"
      )
    ),
    finally = setwd(old_working_directory)
  )
  compile_status <- attr(compile_output, "status")
  if (is.null(compile_status)) {
    compile_status <- 0L
  }
  if (compile_status != 0L || !file.exists(shared_library)) {
    stop(
      paste(
        "Compiling the standalone cloglike source failed:",
        paste(compile_output, collapse = "\n"),
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  set.seed(20260819)
  n <- 80L
  z <- seq(-1.2, 1.2, length.out = n)
  eta <- 0.4 + 0.9 * z
  kappa <- 12
  base_angle <- rnorm(n, sd = 1 / sqrt(kappa))
  y <- 2 * atan(tan(base_angle / 2) + eta)
  Y <- INLA::inla.mdata(y)

  cloglike <- INLA::inla.cloglike.define(
    model = "INLAcirc_cloglike_lavm",
    shlib = normalizePath(shared_library, mustWork = TRUE),
    lavm.link = 0,
    lavm.prior = 0,
    lavm.u = 0.5,
    lavm.alpha = 0.5,
    lavm.initial.theta = log(kappa),
    lavm.fixed.theta = 1
  )

  fit <- INLA::inla(
    Y ~ z,
    family = "cloglike",
    data = list(Y = Y, z = z),
    control.family = list(cloglike = cloglike),
    control.inla = list(
      cmin = 0,
      compute.initial.values = TRUE
    ),
    verbose = FALSE
  )

  coefficient_names <- c("(Intercept)", "z")
  stopifnot(
    inherits(fit, "inla"),
    all(coefficient_names %in% rownames(fit$summary.fixed)),
    all(is.finite(fit$summary.fixed[coefficient_names, "mean"]))
  )

  message("Standalone src-cloglike R/INLA integration test passed.")
  invisible(fit)
}

.INLAcirc_run_test()
