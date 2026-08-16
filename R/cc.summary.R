#' Internal function to compute transformed summaries
#' @keywords internal
cc.summary <- function(result) {

  get_inla_prior <- function(h_name, res) {
    format_pr <- function(pr) {
      if (is.null(pr) || is.null(pr$prior)) return("None")
      if (grepl("expression", pr$prior, ignore.case = TRUE)) return("custom_prior")

      if (pr$prior == "normal") {
        sd_val <- if (pr$param[2] > 0) 1 / sqrt(pr$param[2]) else Inf
        return(sprintf("N(%g, %g)", pr$param[1], signif(sd_val, 4)))
      }
      return(paste0(pr$prior, "(", paste(signif(pr$param, 4), collapse = ", "), ")"))
    }

    if (grepl("observations|von Mises", h_name, ignore.case = TRUE)) {
      for (fam in res$all.hyper$family) {
        fam_base <- sub(" observations.*", "", fam$label)
        if (grepl(fam_base, h_name, ignore.case = TRUE) || grepl(fam$label, h_name, ignore.case = TRUE)) {
          if (length(fam$hyper) > 0) return(format_pr(fam$hyper[[1]]))
        }
      }
    }

    for (rnd in res$all.hyper$random) {
      if (grepl(paste0(" for ", rnd$label, "$"), h_name) || grepl(paste0(" for ", rnd$label, " observations"), h_name)) {
        if (length(rnd$hyper) > 0) {
          if (grepl("^Beta", h_name, ignore.case = TRUE)) {
            for (hyp in rnd$hyper) if (hyp$name == "beta") return(format_pr(hyp))
          } else if (grepl("^Precision|^SD", h_name, ignore.case = TRUE)) {
            for (hyp in rnd$hyper) {
              name_to_check <- if (!is.null(hyp$name)) hyp$name else if (!is.null(hyp$short.name)) hyp$short.name else ""
              if (grepl("prec|theta", name_to_check, ignore.case = TRUE)) return(format_pr(hyp))
            }
          }
          return(format_pr(rnd$hyper[[1]]))
        }
      }
    }

    if (grepl("Precision|SD", h_name, ignore.case = TRUE)) return("loggamma(1, 5e-05)")
    if (grepl("Beta", h_name, ignore.case = TRUE)) return("N(0, 31.62)")
    return("N(0, 31.62)")
  }

  # --- 1. Fixed Effects Processing ---
  fixed_par <- result$summary.fixed[, c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode")]
  fixed_par$prior <- NA

  fe_map <- result$inlacc_meta$fixef_map
  priors_meta <- result$inlacc_meta$priors
  orig_rn <- rownames(fixed_par)

  if (!is.null(fe_map)) {
    new_rn <- orig_rn
    for (i in seq_along(new_rn)) {
      rn <- orig_rn[i]
      mu <- if (!is.null(priors_meta$fixed_means[[rn]])) priors_meta$fixed_means[[rn]] else priors_meta$default_mean
      sd_val <- if (!is.null(priors_meta$fixed_sds[[rn]])) priors_meta$fixed_sds[[rn]] else priors_meta$default_sd
      fixed_par$prior[i] <- sprintf("N(%g, %g)", mu, signif(sd_val, 4))
      if (!is.null(fe_map[[ rn ]])) new_rn[i] <- fe_map[[ rn ]]
    }
    rownames(fixed_par) <- new_rn
  }

  # --- 2. Hyperparameters Processing ---
  hyper_df <- result$summary.hyperpar[, c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode")]
  hyper_df$prior <- NA

  hyper_map <- result$inlacc_meta$hyper_map
  hyper_priors <- result$inlacc_meta$hyper_priors
  lat_processes <- result$inlacc_meta$latent_processes
  meta <- result$inlacc_meta

  for (h_name in rownames(hyper_df)) {
    if (!is.null(hyper_priors) && !is.null(hyper_priors[[h_name]])) {
      hyper_df[h_name, "prior"] <- hyper_priors[[h_name]]
    } else {
      hyper_df[h_name, "prior"] <- get_inla_prior(h_name, result)
    }
  }

  obs_queue <- list()
  for (i in seq_along(meta$responses)) {
    fam <- meta$base_families[i]
    resp <- meta$responses[i]

    if (fam == "lavm") {
      # Outputs \u03BA_ instead of log(\u03BA_
      obs_queue[[length(obs_queue) + 1]] <- list(fam = "lavm", tgt = paste0("\u03BA_", resp), type = "kappa")
    } else if (fam == "gaussian") {
      obs_queue[[length(obs_queue) + 1]] <- list(fam = "gaussian", tgt = paste0("SD_", resp), type = "sd")
    } else if (fam == "gamma") {
      obs_queue[[length(obs_queue) + 1]] <- list(fam = "gamma", tgt = paste0("Shape_", resp), type = "native")
    } else {
      obs_queue[[length(obs_queue) + 1]] <- list(fam = fam, tgt = paste0("Precision_", resp), type = "native")
    }
  }

  new_hyper_lik <- list()
  new_hyper_lat <- list()

  for (h_name in rownames(hyper_df)) {

    base_var <- sub("^.* for ", "", h_name)
    proc_val <- if (!is.null(lat_processes[[base_var]])) lat_processes[[base_var]] else NA

    if (!is.null(hyper_map) && !is.null(hyper_map[[h_name]])) {
      tgt_name <- hyper_map[[h_name]]
      new_row <- hyper_df[h_name, , drop = FALSE]
      new_row$type <- proc_val
      rownames(new_row) <- tgt_name
      new_hyper_lat[[tgt_name]] <- new_row
      next
    }

    # Catch the newly renamed kappa
    if (grepl("observations|von Mises|Theta[0-9]+ for (INLA\\.Data|cloglike)|kappa for lavm", h_name, ignore.case = TRUE)) {
      match_idx <- -1
      if (length(obs_queue) > 0) {
        for (k in seq_along(obs_queue)) {
          expected_fam <- obs_queue[[k]]$fam
          regex_fam <- if (expected_fam == "lavm") "(lavm|von Mises|cloglike|INLA\\.Data)" else expected_fam
          if (grepl(regex_fam, h_name, ignore.case = TRUE)) { match_idx <- k; break }
        }
      }

      if (match_idx > 0) {
        tgt_name <- obs_queue[[match_idx]]$tgt
        obs_type <- obs_queue[[match_idx]]$type
        obs_fam <- obs_queue[[match_idx]]$fam
        obs_queue[[match_idx]] <- NULL

        if (obs_type == "sd") {
          internal_name <- sub("^(Precision|Precision parameter) for ", "Log precision for ", h_name)
          if (internal_name %in% names(result$internal.marginals.hyperpar)) {
            marg <- result$internal.marginals.hyperpar[[internal_name]]
            marg_sd <- INLA::inla.tmarginal(function(x) 1/sqrt(exp(x)), marg)
            z <- INLA::inla.zmarginal(marg_sd, silent = TRUE)
            m <- INLA::inla.mmarginal(marg_sd)
            row_sd <- data.frame(mean = z$mean, sd = z$sd, `0.025quant` = z$quant0.025, `0.5quant` = z$quant0.5, `0.975quant` = z$quant0.975, mode = m, check.names = FALSE)
            row_sd$prior <- hyper_df[h_name, "prior"]
            row_sd$family <- obs_fam
            rownames(row_sd) <- tgt_name
            new_hyper_lik[[tgt_name]] <- row_sd
          } else {
            new_row <- hyper_df[h_name, , drop = FALSE]
            new_row$family <- obs_fam
            rownames(new_row) <- tgt_name
            new_hyper_lik[[tgt_name]] <- new_row
          }
        } else {
          # Handles both 'native' and 'kappa', since lavm is already exponentiated by lavm.rename.inla.output
          new_row <- hyper_df[h_name, , drop = FALSE]
          new_row$family <- obs_fam
          rownames(new_row) <- tgt_name
          new_hyper_lik[[tgt_name]] <- new_row
        }
      } else {
        clean_name <- sub(" observations(\\[[0-9]+\\])?$", "", h_name, ignore.case = TRUE)
        new_row <- hyper_df[h_name, , drop = FALSE]
        new_row$family <- NA
        rownames(new_row) <- clean_name
        new_hyper_lik[[clean_name]] <- new_row
      }

    } else if (grepl("^(Precision|Precision parameter) for ", h_name)) {
      internal_name <- sub("^(Precision|Precision parameter) for ", "Log precision for ", h_name)
      if (internal_name %in% names(result$internal.marginals.hyperpar)) {
        marg <- result$internal.marginals.hyperpar[[internal_name]]
        marg_sd <- INLA::inla.tmarginal(function(x) 1/sqrt(exp(x)), marg)
        z <- INLA::inla.zmarginal(marg_sd, silent = TRUE)
        m <- INLA::inla.mmarginal(marg_sd)
        row_sd <- data.frame(mean = z$mean, sd = z$sd, `0.025quant` = z$quant0.025, `0.5quant` = z$quant0.5, `0.975quant` = z$quant0.975, mode = m, check.names = FALSE)
        row_sd$prior <- hyper_df[h_name, "prior"]
        row_sd$type <- proc_val
        tgt_name <- sub("^(Precision|Precision parameter) for ", "SD_", h_name)
        rownames(row_sd) <- tgt_name
        new_hyper_lat[[tgt_name]] <- row_sd
      } else {
        new_row <- hyper_df[h_name, , drop = FALSE]
        new_row$type <- proc_val
        new_hyper_lat[[h_name]] <- new_row
      }

    } else {
      clean_name <- sub("^Beta(_intern)? for ", "Beta_", h_name)
      clean_name <- sub("^PACF([0-9]+) for ", "PACF\\1_", clean_name)

      new_row <- hyper_df[h_name, , drop = FALSE]
      new_row$type <- proc_val
      rownames(new_row) <- clean_name
      new_hyper_lat[[clean_name]] <- new_row
    }
  }

  lik_out <- do.call(rbind, new_hyper_lik)
  lat_out <- do.call(rbind, new_hyper_lat)

  rename_cols <- function(df) {
    if (is.null(df)) return(NULL)
    cnames <- colnames(df)
    cnames[cnames == "sd"] <- "SD"
    cnames[cnames == "0.025quant"] <- "2.5%"
    cnames[cnames == "0.5quant"] <- "50%"
    cnames[cnames == "0.975quant"] <- "97.5%"
    colnames(df) <- cnames
    return(df)
  }

  fixed_par <- rename_cols(fixed_par)
  lik_out <- rename_cols(lik_out)
  lat_out <- rename_cols(lat_out)

  # Extract CPU tracking
  time_sum <- result$cpu.used
  class(time_sum) <- c("inlacc_time", class(time_sum))

  # Extract Metrics
  metrics <- list()
  if (!is.null(result$mlik)) metrics$mlik <- result$mlik[1,1]

  if (!is.null(result$dic) && !is.null(result$dic$dic)) {
    metrics$dic <- list(
      dic = result$dic$dic,
      peff = result$dic$p.eff
    )
  }
  if (!is.null(result$waic) && !is.null(result$waic$waic)) {
    metrics$waic <- list(
      waic = result$waic$waic,
      peff = result$waic$p.eff
    )
  }

  needs_pred_warning <- FALSE
  if (!is.null(result$inlacc_meta)) {
    has_gamma <- "gamma" %in% result$inlacc_meta$base_families
    has_copied_cov <- length(result$inlacc_meta$copied_covariates) > 0
    has_compute <- isTRUE(result$inlacc_meta$compute_predictor)

    if ((has_gamma || has_copied_cov) && !has_compute) {
      needs_pred_warning <- TRUE
    }
  }

  out <- list(fixed_par = fixed_par,
              likelihood_par = lik_out,
              random_par = lat_out,
              time.summary = time_sum,
              metrics = metrics,
              needs_pred_warning = needs_pred_warning)

  class(out) <- "inlacc_summary"
  return(out)
}

#' Print method for inlacc time summary
#'
#' @param x An `inlacc_time` object.
#' @param ... Additional arguments, currently unused.
#' @return `x`, invisibly.
#' @export
print.inlacc_time <- function(x, ...) {
  cat(sprintf("Setup: %.2f  |  Computation: %.2f  |  Post-processing: %.2f  |  Total: %.2f\n",
              x["Pre"], x["Running"], x["Post"], x["Total"]))
  invisible(x)
}

#' Summary method for inlacc models
#'
#' @param object A fitted `inlacc` model.
#' @param decimal Number of decimal places retained for printing.
#' @param ... Additional arguments, currently unused.
#' @return An object of class `inlacc_summary`.
#' @export
summary.inlacc <- function(object, decimal = 3L, ...) {
  # Generate the base summary using your internal function
  res <- cc.summary(object)

  # Store the user's decimal preference in the object
  res$decimal <- decimal
  return(res)
}

#' Print method for inlacc summaries
#'
#' @param x An `inlacc_summary` object.
#' @param decimal Optional number of decimal places to print.
#' @param metrics Logical; whether to print available model metrics.
#' @param ... Additional arguments, currently unused.
#' @return `x`, invisibly.
#' @export
print.inlacc_summary <- function(x, decimal = NULL, metrics = TRUE, ...) {

  # If decimal isn't provided directly to print(), pull it from the summary object, defaulting to 3
  if (is.null(decimal)) {
    decimal <- if (!is.null(x$decimal)) x$decimal else 3L
  }

  cat("========================================================================================\n")
  cat(" INLAcircular: Bayesian Joint Circular Regression Summary\n")
  cat("========================================================================================\n\n")

  time_sum <- x$time.summary
  if (!is.null(time_sum)) {
    cat("Time Summary (seconds):\n  ")
    print(time_sum)
    cat("\n")
  }

  standardize_df <- function(df, section_name) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$Parameter <- rownames(df)

    if ("family" %in% colnames(df)) {
      df$extra <- as.character(df$family)
    } else if ("type" %in% colnames(df)) {
      df$extra <- as.character(df$type)
    } else {
      df$extra <- ""
    }

    df$Section <- section_name
    df <- df[, c("Section", "Parameter", "mean", "SD", "2.5%", "50%", "97.5%", "mode", "prior", "extra")]
    return(df)
  }

  df_lik <- standardize_df(x$likelihood_par, "Likelihood Parameters")
  df_fix <- standardize_df(x$fixed_par, "Fixed Effect Parameters")
  df_lat <- standardize_df(x$random_par, "Random Effect Parameters")

  df_all <- rbind(df_lik, df_fix, df_lat)
  if (is.null(df_all) || nrow(df_all) == 0) {
    cat("No parameters found.\n")
    return(invisible(x))
  }

  num_cols <- c("mean", "SD", "2.5%", "50%", "97.5%", "mode")
  for (col in num_cols) {
    df_all[[col]] <- format(round(df_all[[col]], decimal), nsmall = decimal, trim = TRUE)
  }

  df_all$prior <- as.character(df_all$prior)
  df_all$prior[is.na(df_all$prior)] <- ""
  df_all$extra[is.na(df_all$extra)] <- ""
  df_all[["|"]] <- "|"

  has_extra <- any(df_all$extra != "")
  print_cols <- c("Parameter", num_cols, "|", "prior")
  if (has_extra) print_cols <- c(print_cols, "extra")

  formatted_headers <- character(length(print_cols))
  for (j in seq_along(print_cols)) {
    col <- print_cols[j]

    if (col == "extra") {
      combined <- c("", "family", "type", df_all[[col]])
      combined_fmt <- format(combined, justify = "left")
      formatted_headers[j] <- combined_fmt[1]
      df_all[[col]] <- combined_fmt[-c(1, 2, 3)]
    } else {
      combined <- c(col, df_all[[col]])
      if (col %in% c("Parameter", "prior")) {
        combined_fmt <- format(combined, justify = "left")
      } else if (col == "|") {
        combined_fmt <- format(combined, justify = "centre")
      } else {
        combined_fmt <- format(combined, justify = "right")
      }
      formatted_headers[j] <- combined_fmt[1]
      df_all[[col]] <- combined_fmt[-1]
    }
  }

  header_line <- paste(formatted_headers, collapse = "  ")

  n_dashes <- nchar(header_line, type = "width")
  if (is.na(n_dashes) || n_dashes < 10) n_dashes <- nchar(header_line)
  sep_line <- strrep("-", n_dashes)

  cat("Posterior Summary :\n")
  cat(sep_line, "\n")
  cat(header_line, "\n")
  cat(sep_line, "\n")

  left_w <- n_dashes
  if (has_extra) {
    left_part <- paste(formatted_headers[print_cols != "extra"], collapse = "  ")
    left_w <- nchar(left_part, type = "width")
    if (is.na(left_w)) left_w <- nchar(left_part)
    left_w <- left_w + 2
  }

  sections <- unique(df_all$Section)
  for (i in seq_along(sections)) {
    sec <- sections[i]

    if (sec == "Likelihood Parameters") {
      sec_title <- "Likelihood Parameters :"
      sec_head <- if (has_extra) "family" else ""
    } else if (sec == "Random Effect Parameters") {
      sec_title <- "Random Effect Parameters :"
      sec_head <- if (has_extra) "type" else ""
    } else {
      sec_title <- paste(sec, ":")
      sec_head <- ""
    }

    if (has_extra && sec_head != "") {
      cat(sprintf("%-*s%s\n", left_w, sec_title, sec_head))
    } else {
      cat(sec_title, "\n")
    }

    sec_rows <- which(df_all$Section == sec)
    for (r in sec_rows) {
      # Print the standard numeric row
      row_vals <- sapply(print_cols, function(col) df_all[[col]][r])
      cat(paste(row_vals, collapse = "  "), "\n")

      # --- INJECT RATE EXPLANATION WITH FAMILY ALIGNMENT ---
      param_clean <- trimws(df_all$Parameter[r])
      if (grepl("^Shape_", param_clean)) {
        resp_name <- sub("^Shape_", "", param_clean)
        msg <- sprintf("[Observation-specific: computed as Shape_%s / exp(linear predictor)]", resp_name)

        param_width <- nchar(df_all$Parameter[r])
        left_string <- sprintf("%-*s  %s", param_width, paste0("Rate_", resp_name), msg)

        if (has_extra) {
          pad_len <- left_w - nchar(left_string, type = "width")
          if (is.na(pad_len)) pad_len <- left_w - nchar(left_string)

          if (pad_len > 0) {
            # Pad exactly to the extra column
            cat(sprintf("%s%s%s\n", left_string, strrep(" ", pad_len), df_all$extra[r]))
          } else {
            # Fallback if the string somehow overflows
            cat(sprintf("%s  %s\n", left_string, df_all$extra[r]))
          }
        } else {
          cat(left_string, "\n")
        }
      }
    }

    cat(sep_line, "\n")
  }

  cat("\n")

  # --- Print Model Performance Metrics ---
  if (metrics && !is.null(x$metrics) && length(x$metrics) > 0) {
    cat("Model Performance Metrics :\n")
    cat(sep_line, "\n")

    if (!is.null(x$metrics$mlik)) {
      cat(sprintf("%-40s: %.*f\n", "Marginal log-Likelihood", decimal, x$metrics$mlik))
    }
    if (!is.null(x$metrics$dic)) {
      cat(sprintf("%-40s: %.*f  (Effective parameters: %.*f)\n",
                  "Deviance Information Criterion (DIC)",
                  decimal, x$metrics$dic$dic, decimal, x$metrics$dic$peff))
    }
    if (!is.null(x$metrics$waic)) {
      cat(sprintf("%-40s: %.*f  (Effective parameters: %.*f)\n",
                  "Watanabe-Akaike Info Criterion (WAIC)",
                  decimal, x$metrics$waic$waic, decimal, x$metrics$waic$peff))
    }
    cat(sep_line, "\n\n")
  }

  # --- INJECT WARNING AT THE VERY END (ONCE) ---
  if (isTRUE(x$needs_pred_warning)) {
    cat("* Note: Add `control.predictor = list(compute = TRUE)` to the inlacc() function\n")
    cat("        to obtain the posterior summary and marginal for the linear predictor.\n\n")
  }

  invisible(x)
}
