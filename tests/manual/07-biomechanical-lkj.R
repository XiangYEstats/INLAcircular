# Developer test: six-response biomechanical model with iidkd_LKJ

required_packages <- c("INLA", "INLAcircular", "circular", "graphpcor")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(INLA))
suppressPackageStartupMessages(library(INLAcircular))

# Locate this file when it is sourced, run with Rscript, or evaluated line by
# line in an open RStudio source editor. No working-directory setup is needed.
.manual_test_directory <- function() {
  is_manual_test_directory <- function(path) {
    is.character(path) &&
      length(path) == 1L &&
      nzchar(path) &&
      dir.exists(file.path(path, "data", "biomechanical"))
  }

  # source("07-biomechanical-lkj.R") or the RStudio Source button.
  for (frame in rev(sys.frames())) {
    path <- frame$ofile
    if (is.character(path) && length(path) == 1L && nzchar(path) &&
        file.exists(path)) {
      candidate <- dirname(normalizePath(path, mustWork = TRUE))
      if (is_manual_test_directory(candidate)) {
        return(candidate)
      }
    }
  }

  # Rscript tests/manual/07-biomechanical-lkj.R.
  file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_argument)) {
    path <- sub("^--file=", "", file_argument[[1L]])
    if (file.exists(path)) {
      candidate <- dirname(normalizePath(path, mustWork = TRUE))
      if (is_manual_test_directory(candidate)) {
        return(candidate)
      }
    }
  }

  # Ctrl+Enter/Run while this file is open in RStudio.
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    context <- tryCatch(
      rstudioapi::getSourceEditorContext(),
      error = function(e) NULL
    )
    path <- if (is.null(context)) NULL else context$path
    if (is.character(path) && length(path) == 1L && nzchar(path) &&
        file.exists(path)) {
      candidate <- dirname(normalizePath(path, mustWork = TRUE))
      if (is_manual_test_directory(candidate)) {
        return(candidate)
      }
    }
  }

  # Common RStudio Project locations, including the package root and its
  # parent directory.
  working_directory <- normalizePath(getwd(), mustWork = TRUE)
  roots <- unique(c(
    working_directory,
    dirname(working_directory),
    dirname(dirname(working_directory))
  ))
  candidates <- unique(c(
    roots,
    file.path(roots, "tests", "manual"),
    file.path(roots, "INLAcircular", "tests", "manual")
  ))
  for (candidate in candidates) {
    if (is_manual_test_directory(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }

  stop(
    paste(
      "Could not locate tests/manual/data/biomechanical.",
      "Open 07-biomechanical-lkj.R in RStudio and run it from that editor,",
      "click Source, or run it with Rscript."
    ),
    call. = FALSE
  )
}

manual_test_dir <- .manual_test_directory()
biomechanical_dir <- file.path(manual_test_dir, "data", "biomechanical")

required_files <- c(
  "rotational_tx_var.csv",
  "rotational_ty_var.csv",
  "rotational_tz_var.csv",
  "translational_x_var.csv",
  "translational_y_var.csv",
  "translational_z_var.csv"
)
missing_files <- required_files[
  !file.exists(file.path(biomechanical_dir, required_files))
]
if (length(missing_files)) {
  stop(
    "Missing biomechanical data file(s): ",
    paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

read_biomechanical <- function(name) {
  utils::read.csv(file.path(biomechanical_dir, name), sep = ";")
}

rot_tx <- read_biomechanical("rotational_tx_var.csv")
rot_ty <- read_biomechanical("rotational_ty_var.csv")
rot_tz <- read_biomechanical("rotational_tz_var.csv")
trans_x <- read_biomechanical("translational_x_var.csv")
trans_y <- read_biomechanical("translational_y_var.csv")
trans_z <- read_biomechanical("translational_z_var.csv")

to_radians <- function(x) {
  as.numeric(circular::conversion.circular(
    circular::circular(x, units = "degrees"),
    units = "radians"
  ))
}

make_configuration <- function(configuration) {
  data.frame(
    run = rot_tx$run,
    load = rot_tx$load,
    loading = ifelse(
      rot_tx$load == "No load",
      0,
      suppressWarnings(as.numeric(gsub("[^0-9]", "", rot_tx$load)))
    ),
    add.load = as.integer(grepl("^L", rot_tx$load)),
    no.load = as.integer(rot_tx$load == "No load"),
    x1 = to_radians(rot_tx[[configuration]]),
    x2 = to_radians(rot_ty[[configuration]]),
    x3 = to_radians(rot_tz[[configuration]]),
    y1 = trans_x[[configuration]],
    y2 = trans_y[[configuration]],
    y3 = trans_z[[configuration]]
  )
}

# This matches configuration c1 in the original biomechanical analysis.
c1 <- rbind(
  make_configuration("OFc1"),
  make_configuration("SNc1")
)
n <- nrow(c1)

lavm_setting <- list(
  link = "inverse.tangent",
  hyper = list(
    kappa = list(
      # LAvM `initial` is already on INLA's internal log(kappa) scale.
      initial = 15,
      fixed = TRUE
    )
  )
)
gaussian_setting <- list(
  hyper = list(
    prec = list(
      initial = 15,
      fixed = TRUE
    )
  )
)

model_c1 <- list(
  likelihood(
    x1 ~ intercept(name = "alpha1") +
      f(i, model = "iidkd_LKJ", pc.prior.u = 1, pc.prior.alpha = 0.5),
    family = "lavm",
    family.setting = lavm_setting
  ),
  likelihood(
    x2 ~ intercept(name = "alpha2") +
      f(i, model = "iidkd_LKJ", pc.prior.u = 1, pc.prior.alpha = 0.5),
    family = "lavm",
    family.setting = lavm_setting
  ),
  likelihood(
    x3 ~ intercept(name = "alpha3") +
      f(i, model = "iidkd_LKJ", pc.prior.u = 1, pc.prior.alpha = 0.5),
    family = "lavm",
    family.setting = lavm_setting
  ),
  likelihood(
    y1 ~ intercept(name = "beta1") +
      f(i, model = "iidkd_LKJ", pc.prior.u = 1, pc.prior.alpha = 0.5),
    family = "gaussian",
    family.setting = gaussian_setting
  ),
  likelihood(
    y2 ~ intercept(name = "beta2") +
      f(i, model = "iidkd_LKJ", pc.prior.u = 1, pc.prior.alpha = 0.5),
    family = "gaussian",
    family.setting = gaussian_setting
  ),
  likelihood(
    y3 ~ intercept(name = "beta3") +
      f(i, model = "iidkd_LKJ", pc.prior.u = 1, pc.prior.alpha = 0.5),
    family = "gaussian",
    family.setting = gaussian_setting
  )
)

fit_c1 <- inlacc(
  model = model_c1,
  data = c1,
  LKJ.eta = 5,
  control.inla = list(compute.initial.values = TRUE),
  control.predictor = list(compute = TRUE),
  metrics = TRUE,
  verbose = TRUE
)

summary(fit_c1)
fit_c1$summary.fixed

# Reconstruct posterior-mode standard deviations and correlations.
theta_i <- fit_c1$lkj_modes[["i"]]
expected_theta <- 6L * (6L + 1L) / 2L
if (length(theta_i) != expected_theta) {
  stop("Could not identify all 21 LKJ hyperparameters in the fitted model.")
}

sds_c1 <- exp(theta_i[seq_len(6L)])
cor_c1 <- graphpcor::basecor(
  base = theta_i[-seq_len(6L)],
  p = 6L
)$base

print(round(sds_c1, 4))
print(round(cor_c1, 2))

# Stacking order: x1, x2, x3, y1, y2, y3.
observed <- as.matrix(c1[c("x1", "x2", "x3", "y1", "y2", "y3")])
fitted <- matrix(
  fit_c1$summary.fitted.values$mean[seq_len(6L * n)],
  nrow = n,
  ncol = 6L
)
mse <- colMeans((observed - fitted)^2)
names(mse) <- colnames(observed)
print(mse)
