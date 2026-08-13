#' Joint Circular Regression with INLA
#'
#' @param model A list of 'likelihood' objects, or a single 'likelihood' object.
#' @param data A data.frame containing all covariates and responses.
#' @param latent.index Optional list of index() mappings.
#' @param LKJ.eta Positive scalar LKJ shape parameter used by every
#'   `iidkd_LKJ` component. Defaults to `5`. Within each corresponding `f()`
#'   term, `pc.prior.u` and `pc.prior.alpha` default to `1` and `0.5`,
#'   respectively.
#' @param control.family Optional global list of control.family configurations.
#' @param control.fixed Optional list for fixed effects priors.
#' @param control.inla Optional control list passed to inla().
#' @param metrics Logical. Compute WAIC, DIC, CPO, etc.?
#' @param verbose Logical. Run INLA in verbose mode?
#' @param ... Additional arguments passed to inla().
#' @return A fitted INLA model.
#' @author Xiang Ye \email{xiang.ye@kaust.edu.sa}
#' @export
inlacc <- function(model, data, latent.index = NULL, LKJ.eta = 5,
                   control.family = NULL,
                   control.fixed = NULL,
                   control.inla = list(cmin = 0, compute.initial.values = TRUE),
                   metrics = FALSE,
                   verbose = FALSE, ...) {

  if (inherits(model, "inlacc_model")) model <- list(model)
  if (!is.list(model) || is.null(model[[1]]$formula) || is.null(model[[1]]$family)) {
    stop("The 'model' argument must be a list of 'likelihood()' objects.")
  }

  dots <- list(...)
  compute_pred <- if (!is.null(dots$control.predictor) && isTRUE(dots$control.predictor$compute)) TRUE else FALSE

  formulas <- lapply(model, function(m) m$formula)
  families <- sapply(model, function(m) m$family)

  n <- nrow(data)
  num_models <- length(formulas)
  responses <- sapply(formulas, function(f) as.character(f)[2])

  # =========================================================================
  # --- LKJ COVARIANCE PRE-PROCESSING ---
  # =========================================================================
  lkj_info <- LKJcc(
    formulas = formulas,
    LKJ.eta = LKJ.eta,
    n.obs = n,
    latent.index = latent.index
  )
  formulas <- lkj_info$formulas
  latent.index <- lkj_info$latent.index

  # =========================================================================
  # --- LATENT INDEX INJECTION ---
  # =========================================================================
  hidden_values_list <- list()
  likelihood_indices <- character()

  if (!is.null(latent.index)) {
    if (!is.null(latent.index$var)) latent.index <- list(latent.index)
  }

  all_f_vars <- character()
  extract_f_vars <- function(expr) {
    if (is.call(expr)) {
      if (identical(expr[[1]], quote(f))) {
        all_f_vars <<- c(all_f_vars, as.character(expr[[2]]))
        if ("replicate" %in% names(expr)) all_f_vars <<- c(all_f_vars, as.character(expr[["replicate"]]))
        if ("group" %in% names(expr)) all_f_vars <<- c(all_f_vars, as.character(expr[["group"]]))
      } else {
        for (j in seq_along(expr)) extract_f_vars(expr[[j]])
      }
    }
  }
  for (f_i in formulas) extract_f_vars(f_i)
  all_f_vars <- unique(all_f_vars)

  modify_f_values <- function(expr, var_name, val_name) {
    if (is.atomic(expr) || is.name(expr)) return(expr)
    if (is.call(expr)) {
      if (identical(expr[[1]], quote(f)) && identical(expr[[2]], as.name(var_name))) {
        expr$values <- as.name(val_name)
        return(expr)
      } else {
        for (i in seq_along(expr)) expr[[i]] <- modify_f_values(expr[[i]], var_name, val_name)
        return(expr)
      }
    }
    return(expr)
  }

  if (!is.null(latent.index)) {
    for (l_idx in latent.index) {
      v_name <- l_idx$var
      if (!(v_name %in% all_f_vars)) {
        warning(sprintf("Latent variable '%s' in latent.index is not used in an f() term. Ignoring.", v_name))
        next
      }

      if (is.character(l_idx$data.id) && length(l_idx$data.id) == 1 && l_idx$data.id == "likelihood") {
        likelihood_indices <- c(likelihood_indices, v_name)
      } else {
        data[[v_name]] <- l_idx$data.id
      }

      if (!is.null(l_idx$process.id)) {
        hidden_val_name <- paste0(".vals_", v_name)
        hidden_values_list[[hidden_val_name]] <- l_idx$process.id
        for (i in seq_along(formulas)) {
          formulas[[i]] <- modify_f_values(formulas[[i]], v_name, hidden_val_name)
          model[[i]]$formula <- formulas[[i]]
        }
      }
    }
  }

  missing_indices <- setdiff(
    all_f_vars,
    c(names(data), likelihood_indices, lkj_info$data.vars)
  )
  for (v in missing_indices) data[[v]] <- 1:n

  # =========================================================================
  # --- COPIED COVARIATES PRE-SCAN ---
  # =========================================================================
  eval_env <- list(
    intercept = function(name = NULL, mean = NULL, sd = NULL) {
      list(name = name, mean = mean, sd = sd)
    },
    covariate = function(var, name = NULL, mean = NULL, sd = NULL, initial = NULL, fixed = FALSE, predictor = FALSE) {
      mc <- match.call()
      list(var = deparse(mc[[2]]), name = name, mean = mean, sd = sd, initial = initial, fixed = fixed, predictor = predictor)
    }
  )

  circ_vars <- responses[families == "lavm"]
  all_rhs_vars <- unique(unlist(lapply(formulas, function(f) all.vars(f[[3]]))))
  circ_covariates <- intersect(circ_vars, all_rhs_vars)

  explicit_pred_vars <- character()
  for (f_i in formulas) {
    trms <- labels(terms(f_i))
    cov_trms <- grep("^covariate\\(", trms, value = TRUE)
    for (ct in cov_trms) {
      parsed <- tryCatch(eval(parse(text = ct), envir = eval_env), error = function(e) NULL)
      if (!is.null(parsed) && isTRUE(parsed$predictor)) {
        explicit_pred_vars <- c(explicit_pred_vars, parsed$var)
      }
    }
  }

  explicit_pred_responses <- intersect(explicit_pred_vars, responses)
  copied_covariates <- unique(c(circ_covariates, explicit_pred_responses))

  # =========================================================================
  # --- DATA BLOCK & LIST-BASED RESPONSE CONSTRUCTION ---
  # =========================================================================
  total_blocks <- num_models + length(copied_covariates)
  Y_list <- lapply(1:total_blocks, function(x) as.numeric(rep(NA, n * total_blocks)))

  data_blocks <- list()
  joint_terms <- character()
  current_dummy_idx <- num_models + 1
  dummy_map <- list()

  for(c_var in copied_covariates) {
    dummy_map[[c_var]] <- current_dummy_idx
    current_dummy_idx <- current_dummy_idx + 1
  }

  local_means <- list()
  local_precs <- list()
  local_sds <- list()
  fixef_map <- list()
  hyper_map <- list()
  hyper_priors <- list()
  latent_processes <- list()

  for (i in 1:num_models) {
    resp <- responses[i]
    block_df <- data.frame(.idx_internal = 1:n)

    for (bv in likelihood_indices) block_df[[bv]] <- i

    trms <- labels(terms(formulas[[i]]))
    has_base_intercept <- attr(terms(formulas[[i]]), "intercept") == 1
    has_custom_intercept <- any(grepl("^intercept\\(", trms))
    model_terms <- character()

    if (has_base_intercept && !has_custom_intercept) {
      int_name <- paste0("Intercept_", resp)
      block_df[[int_name]] <- 1
      model_terms <- c(model_terms, int_name)
      fixef_map[[int_name]] <- int_name
    }

    f_trms <- grep("^f\\(", trms, value = TRUE)
    f_vars <- unique(unlist(lapply(f_trms, function(x) all.vars(parse(text=x)))))
    valid_f_vars <- intersect(f_vars, names(data))
    for(v in valid_f_vars) block_df[[v]] <- data[[v]]

    # LKJcc() owns the multivariate component and replicate indexing. Merge it
    # after ordinary f() variables so the generic data copy cannot overwrite it.
    if (length(lkj_info$block.data[[i]]) > 0L) {
      for (nm in names(lkj_info$block.data[[i]])) {
        block_df[[nm]] <- lkj_info$block.data[[i]][[nm]]
      }
    }

    for (term in trms) {
      if (grepl("^f\\(", term)) {
        model_terms <- c(model_terms, term)
        f_extractor <- function(var, model = "iid", copy = NULL, order = NULL, hyper = NULL, ...) {
          model_expr <- substitute(model)
          if (is.character(model_expr)) {
            m_class <- model_expr[1]
          } else if (is.name(model_expr) &&
                     as.character(model_expr) %in% names(lkj_info$models)) {
            m_class <- "iidkd_LKJ"
          } else {
            m_class <- "custom"
          }
          list(var = deparse(substitute(var)), model = m_class, copy = copy, order = order, hyper = hyper)
        }
        extracted_info <- tryCatch(eval(parse(text = term), envir = list(f = f_extractor)), error = function(e) NULL)
        if (!is.null(extracted_info) && is.list(extracted_info)) {
          f_var <- extracted_info$var
          if (!is.null(extracted_info$copy)) {
            latent_processes[[f_var]] <- "scaling"
          } else if (tolower(extracted_info$model) == "ar") {
            latent_processes[[f_var]] <- sprintf("AR(%d)", if(is.null(extracted_info$order)) 1 else extracted_info$order)
          } else {
            latent_processes[[f_var]] <- toupper(extracted_info$model)
          }
          if (!is.null(extracted_info$hyper) && is.list(extracted_info$hyper)) {
            for (h_key in names(extracted_info$hyper)) {
              h_list <- extracted_info$hyper[[h_key]]
              if (is.list(h_list) && !is.null(h_list$prior)) {
                if (h_list$prior == "normal") {
                  sd_val <- if (!is.null(h_list$param) && length(h_list$param) >= 2 && h_list$param[2] > 0) 1/sqrt(h_list$param[2]) else Inf
                  pr_str <- sprintf("N(%g, %g)", h_list$param[1], signif(sd_val, 4))
                } else {
                  pr_str <- paste0(h_list$prior, "(", paste(signif(h_list$param, 4), collapse = ", "), ")")
                }
                if (grepl("prec|theta", h_key, ignore.case = TRUE)) {
                  hyper_priors[[paste0("Precision for ", f_var)]] <- pr_str
                } else if (grepl("beta", h_key, ignore.case = TRUE)) {
                  hyper_priors[[paste0("Beta for ", f_var)]] <- pr_str
                }
              }
            }
          }
        }
      } else if (grepl("^intercept\\(", term)) {
        parsed <- eval(parse(text = term), envir = eval_env)
        iname <- if (is.null(parsed$name)) paste0("Intercept_", resp) else parsed$name
        block_df[[iname]] <- 1
        model_terms <- c(model_terms, iname)
        fixef_map[[iname]] <- iname
        if (!is.null(parsed$mean)) local_means[[iname]] <- parsed$mean
        if (!is.null(parsed$sd)) {
          local_precs[[iname]] <- 1 / (parsed$sd^2)
          local_sds[[iname]] <- parsed$sd
        }
      } else if (grepl("^covariate\\(", term)) {
        parsed <- eval(parse(text = term), envir = eval_env)
        cvar <- parsed$var
        if (cvar %in% copied_covariates) {
          e_name <- paste0("e_", cvar)
          copy_name <- if (is.null(parsed$name)) paste0("beta_", resp, "_", cvar) else parsed$name
          init_val <- if (!is.null(parsed$initial)) parsed$initial else 0
          fixed_val <- if (isTRUE(parsed$fixed)) "TRUE" else "FALSE"
          if (!is.null(parsed$mean) && !is.null(parsed$sd)) {
            param_str <- sprintf('prior="normal", param=c(%f, %f), ', parsed$mean, 1 / (parsed$sd^2))
            hyper_priors[[paste0("Beta for ", copy_name)]] <- sprintf("N(%g, %g)", parsed$mean, parsed$sd)
          } else {
            param_str <- ""
            hyper_priors[[paste0("Beta for ", copy_name)]] <- "N(0, 31.62)"
          }
          copy_formula <- sprintf(
            'f(%s, copy="%s", hyper=list(beta=list(initial=%f, %sfixed=%s)))',
            copy_name, e_name, init_val, param_str, fixed_val
          )
          model_terms <- c(model_terms, copy_formula)
          block_df[[copy_name]] <- 1:n
          display_name <- if (is.null(parsed$name)) paste0(resp, " ~ \u03B7_", cvar) else parsed$name
          hyper_map[[paste0("Beta for ", copy_name)]] <- display_name
          latent_processes[[copy_name]] <- "scaling"
        } else {
          clean_cvar <- gsub("[^A-Za-z0-9_.]", "_", cvar)
          if (is.null(parsed$name)) {
            cname <- paste0(clean_cvar, "_", i)
            display_name <- paste0(resp, " ~ ", cvar)
          } else {
            cname <- parsed$name
            display_name <- parsed$name
          }
          block_df[[cname]] <- as.numeric(eval(parse(text = cvar), envir = data))
          model_terms <- c(model_terms, cname)
          fixef_map[[cname]] <- display_name
          if (!is.null(parsed$mean)) local_means[[cname]] <- parsed$mean
          if (!is.null(parsed$sd)) {
            local_precs[[cname]] <- 1 / (parsed$sd^2)
            local_sds[[cname]] <- parsed$sd
          }
        }
      } else if (term %in% copied_covariates) {
        e_name <- paste0("e_", term)
        copy_name <- paste0("beta_", resp, "_", term)
        copy_formula <- sprintf('f(%s, copy="%s", hyper=list(beta=list(fixed=FALSE)))', copy_name, e_name)
        model_terms <- c(model_terms, copy_formula)
        block_df[[copy_name]] <- 1:n
        hyper_map[[paste0("Beta for ", copy_name)]] <- paste0(resp, " ~ \u03B7_", term)
        hyper_priors[[paste0("Beta for ", copy_name)]] <- "N(0, 31.62)"
        latent_processes[[copy_name]] <- "scaling"
      } else {
        clean_term <- gsub("[^A-Za-z0-9_.]", "_", term)
        new_name <- paste0(clean_term, "_", i)
        fixef_map[[new_name]] <- paste0(resp, " ~ ", term)
        block_df[[new_name]] <- as.numeric(eval(parse(text = term), envir = data))
        model_terms <- c(model_terms, new_name)
      }
    }

    joint_terms <- c(joint_terms, model_terms)
    Y_list[[i]][((i - 1) * n + 1):(i * n)] <- data[[resp]]

    if (resp %in% copied_covariates) {
      dummy_idx <- dummy_map[[resp]]
      dummy_df <- block_df
      e_name <- paste0("e_", resp)
      ee_name <- paste0("ee_", resp)
      dummy_df[[e_name]] <- 1:n
      dummy_df[[ee_name]] <- -1

      Y_list[[dummy_idx]][((dummy_idx - 1) * n + 1):(dummy_idx * n)] <- 0

      iid_term <- sprintf('f(%s, %s, model="iid", values=1:%d, hyper=list(prec=list(initial=-12, fixed=TRUE)))', e_name, ee_name, n)
      joint_terms <- c(joint_terms, iid_term)
      data_blocks[[dummy_idx]] <- dummy_df
    }

    data_blocks[[i]] <- block_df
  }

  data_stack <- dplyr::bind_rows(data_blocks)

  # =========================================================================
  # --- SAFE CGENERIC INJECTION (Fixes the subscript out of bounds crash) ---
  # =========================================================================
  formula_env <- new.env(parent = environment(formulas[[1]]))
  for (nm in names(lkj_info$models)) {
    assign(nm, lkj_info$models[[nm]], envir = formula_env)
  }

  joint_formula_str <- paste("Y ~ -1 +", paste(unique(joint_terms), collapse = " + "))
  joint_formula <- as.formula(joint_formula_str, env = formula_env)

  ctrl_fixed <- if (is.null(control.fixed)) list() else control.fixed
  user_mean <- if (!is.null(ctrl_fixed$mean)) ctrl_fixed$mean else 0
  user_prec <- if (!is.null(ctrl_fixed$prec)) ctrl_fixed$prec else 0.001

  if (length(local_means) > 0) {
    mean_list <- list(default = user_mean)
    for (nm in names(local_means)) mean_list[[nm]] <- local_means[[nm]]
    ctrl_fixed$mean <- mean_list
  } else {
    ctrl_fixed$mean <- user_mean
  }

  if (length(local_precs) > 0) {
    prec_list <- list(default = user_prec)
    for (nm in names(local_precs)) prec_list[[nm]] <- local_precs[[nm]]
    ctrl_fixed$prec <- prec_list
  } else {
    ctrl_fixed$prec <- user_prec
  }

  # =========================================================================
  # --- ORIGINAL CLOGLIKE & FAMILY INJECTION ---
  # =========================================================================
  joint_families <- character()
  ctrl_fam <- list()

  for (i in seq_along(model)) {
    fam <- families[i]

    usr_ctrl <- NULL
    if (!is.null(control.family) && length(control.family) >= i) {
      usr_ctrl <- control.family[[i]]
    } else if (!is.null(model[[i]]$family.setting)) {
      usr_ctrl <- model[[i]]$family.setting
    }

    if (fam == "lavm") {
      joint_families <- c(joint_families, "cloglike")

      clog_list <- list(cloglike = INLAcircular:::lavm.cloglike(usr_ctrl))

      if (!is.null(usr_ctrl)) {
        for (nn in names(usr_ctrl)) {
          if (!(nn %in% c("link", "hyper", "cloglike", "lambda", "u", "alpha"))) {
            clog_list[[nn]] <- usr_ctrl[[nn]]
          }
        }
      }
      ctrl_fam <- append(ctrl_fam, list(clog_list))
      Y_list[[i]] <- INLA::inla.mdata(Y_list[[i]])

      usr_u <- 0.5
      usr_alpha <- 0.5
      if (!is.null(usr_ctrl$hyper$kappa$param)) {
        usr_u <- as.numeric(usr_ctrl$hyper$kappa$param[1])
        usr_alpha <- as.numeric(usr_ctrl$hyper$kappa$param[2])
      }
      usr_prior <- "pc.vminf"
      if (!is.null(usr_ctrl$hyper$kappa$prior)) {
        usr_prior <- .canonical_lavm_pc_prior(
          usr_ctrl$hyper$kappa$prior
        )
      }
      prior_label <- sprintf("%s(%g, %g)", usr_prior, usr_u, usr_alpha)
      hyper_priors[[sprintf("Theta1 for INLA.Data%d", i)]] <- prior_label
      hyper_priors[["Theta1 for cloglike"]] <- prior_label
    } else {
      joint_families <- c(joint_families, fam)
      if (!is.null(usr_ctrl)) {
        ctrl_fam <- append(ctrl_fam, list(usr_ctrl))
      } else {
        ctrl_fam <- append(ctrl_fam, list(list()))
      }
    }
  }

  for (c_var in copied_covariates) {
    joint_families <- c(joint_families, "gaussian")
    ctrl_fam <- append(ctrl_fam, list(
      list(hyper = list(prec = list(initial = 12, fixed = TRUE)))
    ))
  }

  final_Y <- if (length(Y_list) == 1) Y_list[[1]] else Y_list

  # Cgeneric models are safely isolated in formula_env, avoiding the crash.
  inla_data <- c(list(Y = final_Y), as.list(as.data.frame(data_stack)), hidden_values_list)

  ctrl_compute <- list(cpo = FALSE, dic = FALSE, waic = FALSE, config = FALSE)
  if (isTRUE(metrics)) {
    ctrl_compute$cpo <- TRUE
    ctrl_compute$dic <- TRUE
    ctrl_compute$waic <- TRUE
    ctrl_compute$config <- TRUE
  }

  inla_args <- list(
    formula = joint_formula,
    family = joint_families,
    data = inla_data,
    control.fixed = ctrl_fixed,
    control.family = ctrl_fam,
    control.inla = control.inla,
    control.compute = ctrl_compute,
    verbose = verbose
  )

  # Smoothly merge dots
  for (nm in names(dots)) {
    if (nm == "control.compute") {
      for (sub_nm in names(dots$control.compute)) {
        inla_args$control.compute[[sub_nm]] <- dots$control.compute[[sub_nm]]
      }
    } else {
      inla_args[[nm]] <- dots[[nm]]
    }
  }

  result <- do.call(INLA::inla, inla_args)

  result$inlacc_meta <- list(
    formula = joint_formula,
    data_stack = data_stack,
    families = joint_families,
    control.family = ctrl_fam,
    responses = responses,
    base_families = families,
    fixef_map = fixef_map,
    hyper_map = hyper_map,
    hyper_priors = hyper_priors,
    latent_processes = latent_processes,
    copied_covariates = copied_covariates,
    compute_predictor = compute_pred,
    priors = list(
      fixed_means = local_means,
      fixed_sds = local_sds,
      default_mean = user_mean,
      default_sd = 1 / sqrt(user_prec)
    )
  )

  # =========================================================================
  # --- EXTRACT LKJ HYPERPARAMETER MODES BY LATENT INDEX NAME ---
  # =========================================================================
  lkj_modes <- list()
  if (!is.null(result$mode$theta) && !is.null(result$internal.summary.hyperpar)) {
    hnames <- rownames(result$internal.summary.hyperpar)
    for (vname in lkj_info$lkj_vars) {
      idx <- which(endsWith(hnames, paste0(" for ", vname)))
      if (length(idx) > 0) {
        lkj_modes[[vname]] <- result$mode$theta[idx]
      }
    }
  }
  result$lkj_modes <- lkj_modes
  result$inlacc_meta$LKJ <- lkj_info$processes

  class(result) <- c("inlacc", class(result))
  result <- tryCatch(INLAcircular:::lavm.rename.inla.output(result, result$inlacc_meta), error = function(e) result)
  return(result)
}
