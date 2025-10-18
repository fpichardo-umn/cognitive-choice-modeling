# Single Fit Diagnostics Analysis
# Extract and analyze diagnostics from a single subject fit

suppressPackageStartupMessages({
  library(here)
  library(posterior)
  library(dplyr)
})

# Source required helpers
source(file.path(here::here(), "scripts", "diagnostics", "helpers", "diagnostic_helpers.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_diagnostics.R"))

#' Analyze diagnostics for a single fit
#' @param fit Single fit object
#' @param subject_id Subject identifier (optional)
#' @param thresholds Threshold configuration (optional)
#' @param include_plots Whether to generate plot data (default: TRUE)
#' @return List with complete diagnostic analysis
analyze_single_fit <- function(fit, subject_id = NULL, thresholds = NULL, include_plots = TRUE) {
  if (is.null(thresholds)) {
    thresholds <- load_diagnostic_thresholds()
  }
  
  # Extract basic info
  if (is.null(subject_id) && !is.null(fit$subject_list)) {
    subject_id <- fit$subject_list[1]
  }
  
  analysis <- list(
    metadata = extract_fit_metadata(fit, subject_id),
    convergence = analyze_convergence(fit, thresholds),
    sampling = analyze_sampling_diagnostics(fit, thresholds),
    parameters = analyze_parameter_diagnostics(fit, thresholds),
    summary = create_single_fit_summary(fit, subject_id, thresholds)
  )
  
  # Add plot data if requested
  if (include_plots) {
    analysis$plot_data <- prepare_single_fit_plots(fit)
  }
  
  # Overall assessment
  analysis$overall_status <- analysis$summary$status
  analysis$recommendations <- create_recommendations(analysis$summary, "single")
  
  return(analysis)
}

#' Extract metadata from fit object
#' @param fit Fit object
#' @param subject_id Subject ID
#' @return List with metadata
extract_fit_metadata <- function(fit, subject_id) {
  metadata <- list(
    subject_id = subject_id,
    model_name = if (!is.null(fit$model_name)) fit$model_name else "unknown",
    task = if (!is.null(fit$task)) fit$task else "unknown",
    n_chains = if (!is.null(fit$n_chains)) fit$n_chains else dim(fit$draws)[2],
    n_warmup = if (!is.null(fit$n_warmup)) fit$n_warmup else NA,
    n_iter = if (!is.null(fit$n_iter)) fit$n_iter else NA,
    total_samples = if (!is.null(fit$n_iter) && !is.null(fit$n_chains)) {
      fit$n_iter * fit$n_chains
    } else {
      dim(fit$draws)[1] * dim(fit$draws)[2]
    },
    n_params = if (!is.null(fit$n_params)) fit$n_params else dim(fit$draws)[3],
    adapt_delta = if (!is.null(fit$adapt_delta)) fit$adapt_delta else NA,
    max_treedepth = if (!is.null(fit$max_treedepth)) fit$max_treedepth else NA
  )
  
  # Add sampling info if available
  if (!is.null(fit$data_filt_info)) {
    metadata$data_filt <- fit$data_filt_info
  }
  
  return(metadata)
}

#' Analyze convergence diagnostics
#' @param fit Fit object
#' @param thresholds Threshold configuration
#' @return List with convergence analysis
analyze_convergence <- function(fit, thresholds) {
  convergence <- list()
  
  # Get parameter names
  if (!is.null(fit$params)) {
    params <- fit$params
  } else if (!is.null(fit$all_params)) {
    params <- fit$all_params
  } else {
    params <- dimnames(fit$draws)[[3]]
  }
  
  # R-hat analysis
  if (!is.null(fit$diagnostics) && "rhat" %in% colnames(fit$diagnostics)) {
    rhat_values <- fit$diagnostics[, "rhat"]
    # Filter out lp__ if present
    if ("lp__" %in% names(rhat_values)) {
      rhat_values <- rhat_values[names(rhat_values) != "lp__"]
    }
    
    convergence$rhat <- list(
      values = rhat_values,
      min = min(rhat_values, na.rm = TRUE),
      max = max(rhat_values, na.rm = TRUE),
      median = median(rhat_values, na.rm = TRUE),
      mean = mean(rhat_values, na.rm = TRUE),
      n_problematic = sum(rhat_values > thresholds$thresholds$rhat$problematic, na.rm = TRUE),
      n_acceptable = sum(rhat_values > thresholds$thresholds$rhat$acceptable & 
                        rhat_values <= thresholds$thresholds$rhat$problematic, na.rm = TRUE),
      n_excellent = sum(rhat_values <= thresholds$thresholds$rhat$excellent, na.rm = TRUE),
      problematic_params = names(rhat_values)[rhat_values > thresholds$thresholds$rhat$problematic],
      status = if (max(rhat_values, na.rm = TRUE) <= thresholds$thresholds$rhat$excellent) "PASS"
               else if (max(rhat_values, na.rm = TRUE) <= thresholds$thresholds$rhat$acceptable) "WARN"
               else "FAIL"
    )
  }
  
  # ESS analysis
  if (!is.null(fit$diagnostics)) {
    if ("ess_bulk" %in% colnames(fit$diagnostics)) {
      ess_bulk <- fit$diagnostics[, "ess_bulk"]
      # Filter out lp__ if present
      if ("lp__" %in% names(ess_bulk)) {
        ess_bulk <- ess_bulk[names(ess_bulk) != "lp__"]
      }
      # Total post-warmup samples: n_iter (per chain) * n_chains
      total_samples <- fit$n_iter * fit$n_chains
      ess_bulk_ratio <- ess_bulk / total_samples
      
      convergence$ess_bulk <- list(
        values = ess_bulk,
        ratios = ess_bulk_ratio,
        min = min(ess_bulk, na.rm = TRUE),
        max = max(ess_bulk, na.rm = TRUE),
        median = median(ess_bulk, na.rm = TRUE),
        min_ratio = min(ess_bulk_ratio, na.rm = TRUE),
        median_ratio = median(ess_bulk_ratio, na.rm = TRUE),
        # Classify based on absolute ESS, not ratio
        n_critical = sum(ess_bulk < thresholds$thresholds$ess_bulk$critical, na.rm = TRUE),
        n_acceptable = sum(ess_bulk >= thresholds$thresholds$ess_bulk$critical & 
                          ess_bulk < thresholds$thresholds$ess_bulk$acceptable, na.rm = TRUE),
        n_good = sum(ess_bulk >= thresholds$thresholds$ess_bulk$good, na.rm = TRUE),
        problematic_params = names(ess_bulk)[ess_bulk < thresholds$thresholds$ess_bulk$critical],
        # Status based on absolute ESS counts
        status = if (min(ess_bulk, na.rm = TRUE) >= thresholds$thresholds$ess_bulk$good) "PASS"
                 else if (min(ess_bulk, na.rm = TRUE) >= thresholds$thresholds$ess_bulk$acceptable) "WARN"
                 else "FAIL",
        # Add efficiency note
        efficiency_note = if (min(ess_bulk_ratio, na.rm = TRUE) < thresholds$thresholds$ess_ratio$very_inefficient) {
          sprintf("Very inefficient sampling (%.1f%% efficiency)", min(ess_bulk_ratio, na.rm = TRUE) * 100)
        } else NULL
      )
    }
    
    if ("ess_tail" %in% colnames(fit$diagnostics)) {
      ess_tail <- fit$diagnostics[, "ess_tail"]
      # Filter out lp__ if present
      if ("lp__" %in% names(ess_tail)) {
        ess_tail <- ess_tail[names(ess_tail) != "lp__"]
      }
      # Total post-warmup samples: n_iter (per chain) * n_chains
      total_samples <- fit$n_iter * fit$n_chains
      ess_tail_ratio <- ess_tail / total_samples
      
      convergence$ess_tail <- list(
        values = ess_tail,
        ratios = ess_tail_ratio,
        min = min(ess_tail, na.rm = TRUE),
        median = median(ess_tail, na.rm = TRUE),
        min_ratio = min(ess_tail_ratio, na.rm = TRUE),
        median_ratio = median(ess_tail_ratio, na.rm = TRUE),
        # Classify based on absolute ESS, not ratio
        n_critical = sum(ess_tail < thresholds$thresholds$ess_tail$critical, na.rm = TRUE),
        n_acceptable = sum(ess_tail >= thresholds$thresholds$ess_tail$critical & 
                          ess_tail < thresholds$thresholds$ess_tail$acceptable, na.rm = TRUE),
        n_good = sum(ess_tail >= thresholds$thresholds$ess_tail$good, na.rm = TRUE),
        problematic_params = names(ess_tail)[ess_tail < thresholds$thresholds$ess_tail$critical],
        # Status based on absolute ESS counts
        status = if (min(ess_tail, na.rm = TRUE) >= thresholds$thresholds$ess_tail$good) "PASS"
                 else if (min(ess_tail, na.rm = TRUE) >= thresholds$thresholds$ess_tail$acceptable) "WARN"
                 else "FAIL"
      )
    }
  }
  
  # MCSE analysis
  if (!is.null(fit$draws)) {
    draws_matrix <- posterior::as_draws_matrix(fit$draws)
    # Filter out lp__ if present
    if ("lp__" %in% colnames(draws_matrix)) {
      draws_matrix <- draws_matrix[, colnames(draws_matrix) != "lp__", drop = FALSE]
    }
    mcse_values <- apply(draws_matrix, 2, posterior::mcse_mean)
    posterior_sd <- apply(draws_matrix, 2, sd, na.rm = TRUE)
    mcse_ratio <- mcse_values / posterior_sd
    
    convergence$mcse <- list(
      values = mcse_values,
      ratios = mcse_ratio,
      min = min(mcse_values, na.rm = TRUE),
      max = max(mcse_values, na.rm = TRUE),
      median = median(mcse_values, na.rm = TRUE),
      min_ratio = min(mcse_ratio, na.rm = TRUE),
      max_ratio = max(mcse_ratio, na.rm = TRUE),
      median_ratio = median(mcse_ratio, na.rm = TRUE),
      n_problematic = sum(mcse_ratio > thresholds$thresholds$mcse_ratio$problematic, na.rm = TRUE),
      n_acceptable = sum(mcse_ratio > thresholds$thresholds$mcse_ratio$acceptable & 
                        mcse_ratio <= thresholds$thresholds$mcse_ratio$borderline, na.rm = TRUE),
      problematic_params = names(mcse_ratio)[mcse_ratio > thresholds$thresholds$mcse_ratio$problematic],
      status = if (max(mcse_ratio, na.rm = TRUE) <= thresholds$thresholds$mcse_ratio$acceptable) "PASS"
               else if (max(mcse_ratio, na.rm = TRUE) <= thresholds$thresholds$mcse_ratio$borderline) "WARN"
               else "FAIL"
    )
  }
  
  return(convergence)
}

#' Analyze sampling diagnostics
#' @param fit Fit object
#' @param thresholds Threshold configuration
#' @return List with sampling analysis
analyze_sampling_diagnostics <- function(fit, thresholds) {
  sampling <- list()
  
  # Divergences
  if (!is.null(fit$sampler_diagnostics)) {
    divergences <- sum(fit$sampler_diagnostics[,,"divergent__"], na.rm = TRUE)
    n_iter <- dim(fit$sampler_diagnostics)[1]
    n_chains <- dim(fit$sampler_diagnostics)[2]
    total_transitions <- n_iter * n_chains
    divergence_rate <- divergences / total_transitions
    
    sampling$divergences <- list(
      count = divergences,
      rate = divergence_rate,
      total_transitions = total_transitions,
      by_chain = colSums(fit$sampler_diagnostics[,,"divergent__"], na.rm = TRUE),
      status = if (divergence_rate <= thresholds$thresholds$divergence_rate$acceptable) "PASS"
               else if (divergence_rate <= thresholds$thresholds$divergence_rate$borderline) "WARN"
               else "FAIL"
    )
  } else if (!is.null(fit$diagnostic_summary$num_divergent)) {
    divergences <- sum(fit$diagnostic_summary$num_divergent)
    # n_iter is already post-warmup samples per chain
    total_transitions <- fit$n_iter * fit$n_chains
    divergence_rate <- divergences / total_transitions
    
    sampling$divergences <- list(
      count = divergences,
      rate = divergence_rate,
      total_transitions = total_transitions,
      by_chain = fit$diagnostic_summary$num_divergent,
      status = if (divergence_rate <= thresholds$thresholds$divergence_rate$acceptable) "PASS"
               else if (divergence_rate <= thresholds$thresholds$divergence_rate$borderline) "WARN"
               else "FAIL"
    )
  }
  
  # Max treedepth
  if (!is.null(fit$sampler_diagnostics) && "treedepth__" %in% dimnames(fit$sampler_diagnostics)[[3]]) {
    max_td <- if (!is.null(fit$max_treedepth)) fit$max_treedepth else 10
    treedepth_hits <- sum(fit$sampler_diagnostics[,,"treedepth__"] >= max_td, na.rm = TRUE)
    total_transitions <- dim(fit$sampler_diagnostics)[1] * dim(fit$sampler_diagnostics)[2]
    
    sampling$treedepth <- list(
      max_allowed = max_td,
      hits = treedepth_hits,
      rate = treedepth_hits / total_transitions,
      by_chain = colSums(fit$sampler_diagnostics[,,"treedepth__"] >= max_td, na.rm = TRUE),
      status = if (treedepth_hits == 0) "PASS" else if (treedepth_hits < total_transitions * 0.01) "WARN" else "FAIL"
    )
  } else if (!is.null(fit$diagnostic_summary$num_max_treedepth)) {
    treedepth_hits <- sum(fit$diagnostic_summary$num_max_treedepth)
    # n_iter is already post-warmup samples per chain
    total_transitions <- fit$n_iter * fit$n_chains
    
    sampling$treedepth <- list(
      max_allowed = fit$max_treedepth,
      hits = treedepth_hits,
      rate = treedepth_hits / total_transitions,
      by_chain = fit$diagnostic_summary$num_max_treedepth,
      status = if (treedepth_hits == 0) "PASS" else if (treedepth_hits < total_transitions * 0.01) "WARN" else "FAIL"
    )
  }
  
  # EBFMI
  if (!is.null(fit$diagnostic_summary$ebfmi)) {
    ebfmi_values <- fit$diagnostic_summary$ebfmi
    min_ebfmi <- min(ebfmi_values, na.rm = TRUE)
    
    sampling$ebfmi <- list(
      values = ebfmi_values,
      min = min_ebfmi,
      median = median(ebfmi_values, na.rm = TRUE),
      status = if (min_ebfmi >= thresholds$thresholds$ebfmi$acceptable) "PASS" else "FAIL"
    )
  } else if (!is.null(fit$sampler_diagnostics) && "energy__" %in% dimnames(fit$sampler_diagnostics)[[3]]) {
    # Calculate EBFMI manually
    n_chains <- dim(fit$sampler_diagnostics)[2]
    ebfmi_values <- sapply(1:n_chains, function(chain) {
      energy <- fit$sampler_diagnostics[, chain, "energy__"]
      sum(diff(energy)^2) / (length(energy) - 1) / var(energy)
    })
    min_ebfmi <- min(ebfmi_values, na.rm = TRUE)
    
    sampling$ebfmi <- list(
      values = ebfmi_values,
      min = min_ebfmi,
      median = median(ebfmi_values, na.rm = TRUE),
      status = if (min_ebfmi >= thresholds$thresholds$ebfmi$acceptable) "PASS" else "FAIL"
    )
  }
  
  # Energy distribution (for normality check)
  if (!is.null(fit$sampler_diagnostics) && "energy__" %in% dimnames(fit$sampler_diagnostics)[[3]]) {
    energy <- as.vector(fit$sampler_diagnostics[,, "energy__"])
    
    # Test normality if enough samples
    if (length(energy) >= 3) {
      sampling$energy <- list(
        mean = mean(energy, na.rm = TRUE),
        sd = sd(energy, na.rm = TRUE),
        median = median(energy, na.rm = TRUE),
        n_samples = length(energy)
      )
    }
  }
  
  return(sampling)
}

#' Analyze parameter-specific diagnostics
#' @param fit Fit object
#' @param thresholds Threshold configuration
#' @return Data frame with per-parameter diagnostics
analyze_parameter_diagnostics <- function(fit, thresholds) {
  # Get parameter names (exclude lp__)
  if (!is.null(fit$params)) {
    params <- fit$params
  } else if (!is.null(fit$all_params)) {
    params <- fit$all_params
  } else {
    params <- dimnames(fit$draws)[[3]]
  }
  
  # Filter out lp__ if present
  params <- params[params != "lp__"]
  n_params <- length(params)
  
  # Initialize data frame
  param_df <- data.frame(
    parameter = params,
    rhat = rep(NA_real_, n_params),
    rhat_status = rep(NA_character_, n_params),
    ess_bulk = rep(NA_real_, n_params),
    ess_bulk_ratio = rep(NA_real_, n_params),
    ess_bulk_status = rep(NA_character_, n_params),
    ess_tail = rep(NA_real_, n_params),
    ess_tail_ratio = rep(NA_real_, n_params),
    ess_tail_status = rep(NA_character_, n_params),
    mcse = rep(NA_real_, n_params),
    mcse_ratio = rep(NA_real_, n_params),
    mcse_status = rep(NA_character_, n_params),
    stringsAsFactors = FALSE
  )
  
  # Fill in diagnostics (filter out lp__ by matching param names)
  if (!is.null(fit$diagnostics)) {
    if ("rhat" %in% colnames(fit$diagnostics)) {
      rhat_all <- fit$diagnostics[, "rhat"]
      # Match by parameter names, excluding lp__
      rhat_filtered <- rhat_all[names(rhat_all) %in% params]
      param_df$rhat <- rhat_filtered[params]  # Ensure correct order
      param_df$rhat_status <- sapply(param_df$rhat, function(x) {
        classify_diagnostic(x, thresholds, "rhat")
      })
    }
    
    if ("ess_bulk" %in% colnames(fit$diagnostics)) {
      ess_bulk_all <- fit$diagnostics[, "ess_bulk"]
      ess_bulk_filtered <- ess_bulk_all[names(ess_bulk_all) %in% params]
      param_df$ess_bulk <- ess_bulk_filtered[params]
      # Total post-warmup samples: n_iter (per chain) * n_chains
      total_samples <- fit$n_iter * fit$n_chains
      param_df$ess_bulk_ratio <- param_df$ess_bulk / total_samples
      # Use ABSOLUTE ESS for status, not ratio
      param_df$ess_bulk_status <- sapply(param_df$ess_bulk, function(x) {
        classify_diagnostic(x, thresholds, "ess_bulk")
      })
    }
    
    if ("ess_tail" %in% colnames(fit$diagnostics)) {
      ess_tail_all <- fit$diagnostics[, "ess_tail"]
      ess_tail_filtered <- ess_tail_all[names(ess_tail_all) %in% params]
      param_df$ess_tail <- ess_tail_filtered[params]
      # Total post-warmup samples: n_iter (per chain) * n_chains
      total_samples <- fit$n_iter * fit$n_chains
      param_df$ess_tail_ratio <- param_df$ess_tail / total_samples
      # Use ABSOLUTE ESS for status, not ratio
      param_df$ess_tail_status <- sapply(param_df$ess_tail, function(x) {
        classify_diagnostic(x, thresholds, "ess_tail")
      })
    }
  }
  
  # MCSE (also filter lp__ from draws)
  if (!is.null(fit$draws)) {
    draws_matrix <- posterior::as_draws_matrix(fit$draws)
    # Filter out lp__ column
    if ("lp__" %in% colnames(draws_matrix)) {
      draws_matrix <- draws_matrix[, colnames(draws_matrix) != "lp__", drop = FALSE]
    }
    # Match columns to params
    draws_matrix <- draws_matrix[, params, drop = FALSE]
    param_df$mcse <- apply(draws_matrix, 2, posterior::mcse_mean)
    posterior_sd <- apply(draws_matrix, 2, sd, na.rm = TRUE)
    param_df$mcse_ratio <- param_df$mcse / posterior_sd
    param_df$mcse_status <- sapply(param_df$mcse_ratio, function(x) {
      classify_diagnostic(x, thresholds, "mcse_ratio")
    })
  }
  
  return(param_df)
}

#' Prepare plot data for single fit
#' @param fit Fit object
#' @return List with plot data structures
prepare_single_fit_plots <- function(fit) {
  plot_data <- list()
  
  # Traceplot data
  if (!is.null(fit$draws)) {
    # Sample parameters for traceplots (max 10), excluding lp__
    all_params <- dimnames(fit$draws)[[3]]
    all_params <- all_params[all_params != "lp__"]
    n_sample <- min(10, length(all_params))
    sampled_params <- sample(all_params, n_sample)
    
    plot_data$traceplot_params <- sampled_params
  }
  
  # Energy distribution
  if (!is.null(fit$sampler_diagnostics) && "energy__" %in% dimnames(fit$sampler_diagnostics)[[3]]) {
    plot_data$energy <- as.vector(fit$sampler_diagnostics[,, "energy__"])
  }
  
  # Divergence locations (if any)
  if (!is.null(fit$sampler_diagnostics) && "divergent__" %in% dimnames(fit$sampler_diagnostics)[[3]]) {
    divergent <- fit$sampler_diagnostics[,, "divergent__"]
    if (sum(divergent) > 0) {
      plot_data$has_divergences <- TRUE
      plot_data$divergence_iterations <- which(divergent == 1, arr.ind = TRUE)
    } else {
      plot_data$has_divergences <- FALSE
    }
  }
  
  return(plot_data)
}

#' Get problematic parameters
#' @param analysis Analysis results from analyze_single_fit
#' @param thresholds Threshold configuration
#' @return Data frame with problematic parameters only
get_problematic_parameters <- function(analysis, thresholds) {
  param_df <- analysis$parameters
  
  # Filter to problematic parameters
  problematic <- param_df[
    (param_df$rhat_status == "FAIL" | param_df$rhat_status == "WARN") |
    (param_df$ess_bulk_status == "FAIL" | param_df$ess_bulk_status == "WARN") |
    (param_df$mcse_status == "FAIL" | param_df$mcse_status == "WARN"),
  ]
  
  if (nrow(problematic) == 0) {
    return(NULL)
  }
  
  # Sort by severity (FAIL first, then WARN)
  problematic$severity_score <- 
    (problematic$rhat_status == "FAIL") * 10 +
    (problematic$rhat_status == "WARN") * 5 +
    (problematic$ess_bulk_status == "FAIL") * 10 +
    (problematic$ess_bulk_status == "WARN") * 5 +
    (problematic$mcse_status == "FAIL") * 10 +
    (problematic$mcse_status == "WARN") * 5
  
  problematic <- problematic[order(problematic$severity_score, decreasing = TRUE), ]
  problematic$severity_score <- NULL
  
  return(problematic)
}
