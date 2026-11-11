#!/usr/bin/env Rscript

# Render the repository guide using the same function as the editor's Knit button.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  supported <- c("--help", "-h", "--run-models", "--run-lkj")
  unknown <- setdiff(args, supported)
  if (length(unknown)) {
    stop(
      "Unknown argument(s): ", paste(unknown, collapse = ", "),
      ". Run Rscript tools/render-guide.R --help for usage.",
      call. = FALSE
    )
  }
  if (any(args %in% c("--help", "-h"))) {
    cat(
      "Usage: Rscript tools/render-guide.R [--run-models] [--run-lkj]\n\n",
      "Render output/html/INLAcircular-guide.html and\n",
      "output/pdf/INLAcircular-guide.pdf from the current source tree.\n\n",
      "By default, run the lightweight examples and display model-fitting\n",
      "code without executing it, regardless of existing environment flags.\n\n",
      "  --run-models  Also run model examples; requires INLA >= 25.08.21.\n",
      "  --run-lkj     Also run LKJ examples; implies --run-models and\n",
      "                additionally requires graphpcor and INLAtools.\n",
      "  --help, -h    Show this help.\n\n",
      "Requires pkgload, rmarkdown, knitr, Pandoc, and a LaTeX installation\n",
      "providing pdflatex. The script can be called from any directory.\n",
      sep = ""
    )
    return(invisible(NULL))
  }

  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg) != 1L) {
    stop("Run this script with Rscript tools/render-guide.R.", call. = FALSE)
  }
  script <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
  input <- file.path(root, "vignettes", "INLAcircular-guide.Rmd")
  if (!file.exists(input)) {
    stop("Cannot find the guide source: ", input, call. = FALSE)
  }

  run_lkj <- "--run-lkj" %in% args
  run_models <- run_lkj || "--run-models" %in% args
  required <- c("pkgload", "rmarkdown", "knitr")
  if (run_models) required <- c(required, "INLA")
  if (run_lkj) required <- c(required, "graphpcor", "INLAtools")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      ". See the guide's installation instructions before rendering.",
      call. = FALSE
    )
  }
  if (run_models && utils::packageVersion("INLA") < "25.08.21") {
    stop("Model examples require INLA >= 25.08.21.", call. = FALSE)
  }
  if (!rmarkdown::pandoc_available()) {
    stop("Rendering the guide requires Pandoc. Install Pandoc first.", call. = FALSE)
  }
  if (!nzchar(Sys.which("pdflatex"))) {
    stop(
      "Rendering the PDF requires pdflatex. Install a LaTeX distribution ",
      "such as TinyTeX first.",
      call. = FALSE
    )
  }

  flag_names <- c("INLACIRCULAR_RUN_MODELS", "INLACIRCULAR_RUN_LKJ")
  old_flags <- Sys.getenv(flag_names, unset = NA_character_)
  on.exit({
    unset <- is.na(old_flags)
    if (any(unset)) Sys.unsetenv(names(old_flags)[unset])
    if (any(!unset)) do.call(Sys.setenv, as.list(old_flags[!unset]))
  }, add = TRUE)
  Sys.setenv(
    INLACIRCULAR_RUN_MODELS = if (run_models) "true" else "false",
    INLACIRCULAR_RUN_LKJ = if (run_lkj) "true" else "false"
  )
  if (run_models) {
    INLA::inla.setOption(num.threads = "2:1", blas.num.threads = 1L)
  }

  metadata <- rmarkdown::yaml_front_matter(input)
  if (!is.character(metadata$knit) || length(metadata$knit) != 1L) {
    stop("The guide must define its rendering function in YAML 'knit'.", call. = FALSE)
  }
  knit <- eval(parse(text = metadata$knit), envir = globalenv())
  if (!is.function(knit)) {
    stop("The guide's YAML 'knit' entry must evaluate to a function.", call. = FALSE)
  }
  message(
    "Rendering the guide (model examples: ", run_models,
    "; LKJ examples: ", run_lkj, ")."
  )
  invisible(knit(input, encoding = "UTF-8"))
}

main()
