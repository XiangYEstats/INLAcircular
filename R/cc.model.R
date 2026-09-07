#' Define one likelihood block for `inlacc()`
#'
#' A likelihood block stores one response formula, its likelihood family, and
#' the corresponding INLA family controls. Pass one block, or a list of blocks,
#' to [inlacc()].
#'
#' @param formula A two-sided model formula defining one response block. The
#'   left side should be a single column name in `data`; the right side may
#'   contain numeric fixed-effect expressions, [intercept()], [covariate()],
#'   and bare INLA `f()` terms.
#' @param family One character string naming the likelihood family. Use
#'   `"lavm"` for the package's link-adjusted von Mises likelihood. Other
#'   family names are passed to `INLA::inla()`.
#' @param family.setting Optional per-block INLA `control.family` list. For an
#'   LAvM block, use `link` and `hyper$kappa`; see [lavm.cloglike()].
#'
#' @details All blocks supplied to one [inlacc()] call share the rows of a
#' single data frame. Their order determines the stacking order and, for an
#' `iidkd_LKJ` group, the multivariate component order. This constructor stores
#' the specification; validation that depends on the data occurs in
#' [inlacc()].
#'
#' @return An object of class `inlacc_model` containing `formula`, `family`,
#'   and `family.setting`.
#' @seealso [inlacc()], [intercept()], [covariate()], [index()]
#' @examples
#' block <- likelihood(
#'   angle ~ intercept(name = "alpha0", mean = 0, sd = 2) +
#'     covariate(z, name = "alpha_z", mean = 0, sd = 1),
#'   family = "lavm",
#'   family.setting = list(
#'     link = "inverse.tangent",
#'     hyper = list(
#'       kappa = list(
#'         initial = 6,
#'         prior = "pc.vminf",
#'         param = c(0.5, 0.5),
#'         fixed = FALSE
#'       )
#'     )
#'   )
#' )
#' @export
likelihood <- function(formula, family, family.setting = NULL) {
  if (!inherits(formula, "formula")) stop("'formula' must be a valid R formula.")
  if (!is.character(family)) stop("'family' must be a character string (e.g., \"lavm\").")

  structure(list(formula = formula,
                 family = family,
                 family.setting = family.setting),
            class = "inlacc_model")
}

#' Define a custom intercept in an `inlacc()` formula
#'
#' `intercept()` is a formula marker interpreted by [inlacc()]. It replaces
#' the automatically named intercept for that likelihood block.
#'
#' @param name Optional coefficient name. The default is
#'   `Intercept_<response>`. Use unique, syntactically valid names unless a
#'   shared fixed coefficient is intentional.
#' @param mean Optional Gaussian prior mean. If omitted, the corresponding
#'   `control.fixed` mean is inherited.
#' @param sd Optional positive Gaussian prior standard deviation. If omitted,
#'   the corresponding `control.fixed` precision is inherited.
#'
#' @details A normal formula intercept (`1`) creates a block-specific
#' automatic intercept. `0` or `-1` removes it. An explicit `intercept()` term
#' suppresses the automatic intercept and adds the custom term, even in a
#' no-intercept formula. The marker is not intended to be evaluated outside a
#' [likelihood()] formula.
#'
#' @return `NULL`, invisibly; [inlacc()] interprets the unevaluated call.
#' @seealso [likelihood()], [covariate()], [inlacc()]
#' @examples
#' likelihood(
#'   y ~ intercept(name = "beta0", mean = 0, sd = 5),
#'   family = "gaussian"
#' )
#' @export
intercept <- function(name = NULL, mean = NULL, sd = NULL) {
  invisible(NULL)
}

#' Define a covariate or copied predictor in an `inlacc()` formula
#'
#' `covariate()` is a formula marker with two modes. For an ordinary data
#' variable it defines a fixed effect. When `var` is also a response in another
#' block, it can define an INLA copy-scaling coefficient for that response's
#' latent predictor.
#'
#' @param var An unquoted data-column name or numeric expression. This is an
#'   argument name, not a separate `var()` function.
#' @param name Optional coefficient and display name.
#' @param mean Optional Gaussian prior mean.
#' @param sd Optional Gaussian prior standard deviation.
#' @param initial Optional initial value for a copied predictor's INLA scaling
#'   coefficient. It is ignored for an ordinary fixed covariate. The copied
#'   coefficient default is `0`.
#' @param fixed Logical; for a copied predictor, whether to fix the scaling
#'   coefficient at `initial`. It is ignored for an ordinary fixed covariate.
#' @param predictor `NULL` or one logical value. The default `NULL`
#'   automatically uses an INLA copy term when `var` is an untransformed
#'   response from another likelihood block. `TRUE` explicitly requests that
#'   behavior and errors unless `var` is such a response. `FALSE` always uses
#'   the observed value as an ordinary fixed covariate.
#'
#' @details For an ordinary fixed effect, `mean` and `sd` independently
#' override `control.fixed`; omitted values inherit the global settings. For a
#' copied predictor, both must be supplied to customize the scaling prior. If
#' either is missing, the current implementation uses its default
#' `N(0, 31.62^2)` prior. Set
#' `control.predictor = list(compute = TRUE)` in [inlacc()] when fitted latent
#' predictors are required. A response used as a bare formula term is also
#' copied automatically; write `covariate(response, predictor = FALSE)` to use
#' its observed values instead. Transformed response expressions are ordinary
#' observed fixed effects unless `predictor = TRUE`, which is invalid because
#' a transformed latent predictor cannot be copied by this interface.
#'
#' @return `NULL`, invisibly; [inlacc()] interprets the unevaluated call.
#' @seealso [likelihood()], [intercept()], [inlacc()]
#' @examples
#' likelihood(
#'   y ~ covariate(z, name = "beta_z", mean = 0, sd = 2),
#'   family = "gaussian"
#' )
#'
#' # In a joint model, `angle` is copied automatically because it is another
#' # likelihood block's response. Use predictor = TRUE to make that explicit.
#' likelihood(
#'   y ~ covariate(
#'     angle,
#'     name = "beta_angle_predictor",
#'     mean = 0,
#'     sd = 1,
#'     initial = 0,
#'     fixed = FALSE,
#'     predictor = TRUE
#'   ),
#'   family = "gaussian"
#' )
#' @export
covariate <- function(var, name = NULL, mean = NULL, sd = NULL, initial = NULL, fixed = FALSE, predictor = NULL) {
  invisible(NULL)
}

#' Map observations to a latent process in `inlacc()`
#'
#' Create an observation-to-process mapping for an INLA `f()` term. Pass one
#' mapping, or a list of mappings, through the `latent.index` argument of
#' [inlacc()].
#'
#' @param var Character string naming the first argument of the corresponding
#'   bare `f()` term, for example `"hour_effect"`. Quote the name.
#' @param data.id Observation-level node identifiers. This should normally
#'   have `nrow(data)` entries. The special value `"likelihood"` assigns the
#'   likelihood-block number for an ordinary shared effect.
#' @param process.id Optional complete ordered node set. When supplied for a
#'   matching first `f()` index, [inlacc()] injects it as the INLA `values`
#'   argument, replacing an existing `values` argument and allowing unobserved
#'   process nodes. If `NULL`, no explicit `values` vector is injected and
#'   INLA infers the represented values from `data.id`.
#'
#' @details If `var` is not used as an `f()` first, group, or replicate
#' argument, [inlacc()] ignores the mapping with a warning. Without a mapping
#' or same-named data column, a discovered index is generated as
#' `seq_len(nrow(data))`. A mapping may provide `data.id` for a first, group,
#' or replicate variable, but `process.id` is injected only when `var` matches
#' the first argument of `f()`.
#'
#' For `model = "iidkd_LKJ"`, the special form
#' `index(var, data.id = replicate_ids, process.id = "likelihood")` has a
#' different meaning: `data.id` supplies positive integer replicate IDs and
#' the component index is generated from likelihood order. See [LKJcc()].
#'
#' @return A list containing `var`, `data.id`, and `process.id`.
#' @seealso [inlacc()], [LKJcc()]
#' @examples
#' hour_index <- index(
#'   var = "hour_effect",
#'   data.id = c(1, 2, 2, 4),
#'   process.id = 1:4
#' )
#' @export
index <- function(var, data.id, process.id = NULL) {
  list(
    var = as.character(var),
    data.id = data.id,
    process.id = process.id
  )
}
