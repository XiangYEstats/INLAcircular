#' Define a Likelihood Block for INLAcc
#' @export
likelihood <- function(formula, family, family.setting = NULL) {
  if (!inherits(formula, "formula")) stop("'formula' must be a valid R formula.")
  if (!is.character(family)) stop("'family' must be a character string (e.g., \"lavm\").")

  structure(list(formula = formula,
                 family = family,
                 family.setting = family.setting),
            class = "inlacc_model")
}

#' Define Custom Intercept Settings
#' @export
intercept <- function(name = NULL, mean = NULL, sd = NULL) {
  invisible(NULL)
}

#' Define Custom Covariate Settings
#' @export
covariate <- function(var, name = NULL, mean = NULL, sd = NULL, initial = NULL, fixed = FALSE, predictor = FALSE) {
  invisible(NULL)
}

#' Define identifier mapping for latent processes
#'
#' @param var Character string. The name of the latent variable as it appears in your model formula (e.g., "w2").
#' @param data.id A vector mapping each observation in your dataset to a specific location/node in the latent process. The length must match the number of rows in your data.
#' @param process.id A vector defining the complete, unique set of locations/nodes that make up the latent process (e.g., \code{1:24} for hours of a day, or \code{1:50} for spatial regions). If \code{NULL}, it automatically defaults to the sorted unique values found in \code{data.id}.
#' @return A list containing the mapping specifications.
#' @export
index <- function(var, data.id, process.id = NULL) {
  list(
    var = as.character(var),
    data.id = data.id,
    process.id = process.id
  )
}
