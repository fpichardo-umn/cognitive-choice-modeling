# Helper functions for diagnostic analysis
# These functions provide utilities for processing and analyzing MCMC diagnostics

suppressPackageStartupMessages({
  library(dplyr)
  library(posterior)
  library(yaml)
  library(here)
  library(rlang)
})

# Define null-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Format status badge for HTML output
#' @param status Status string ("PASS", "WARN", or "FAIL")
#' @return HTML string for badge
format_status_badge <- function(status) {
  color <- switch(status,
    "PASS" = "#28a745",
    "WARN" = "#ffc107",
    "FAIL" = "#dc3545",
    "#6c757d"  # default gray
  )
  
  sprintf('<span style="display:inline-block;padding:4px 8px;border-radius:4px;background-color:%s;color:white;font-weight:bold;">%s</span>',
          color, status)
}

#' Load diagnostic thresholds from config
#' @param config_file Path to config file (optional)
#' @return List with threshold values
load_diagnostic_thresholds <- function(config_file = NULL) {
  if (is.null(config_file)) {
    config_file <- file.path(here::here(), "scripts", "diagnostics", "config", "diagnostic_thresholds.yaml")
  }
  
  if (!file.exists(config_file)) {
    # Return default thresholds
    return(list(
      thresholds = list(
        rhat = list(excellent = 1.01, acceptable = 1.05, problematic = 1.1),
        ess_ratio = list(good = 0.7, acceptable = 0.5, problematic = 0.3),
        divergence_rate = list(acceptable = 0.001, borderline = 0.01, problematic = 0.05),
        mcse_ratio = list(acceptable = 0.05, borderline = 0.1, problematic = 0.15),
        ebfmi = list(acceptable = 0.3, problematic = 0.2)
      ),
      problem_weights = list(
        high_rhat = 5,
        divergences = 3,
        low_ess = 2,
        high_mcse = 1,
        low_ebfmi = 2
      )
    ))
  }
  
  yaml::read_yaml(config_file)
}

#' Load report configuration
#' @param config_file Path to config file (optional)
#' @return List with report configuration
load_report_config <- function(config_file = NULL) {
  if (is.null(config_file)) {
    config_file <- file.path(here::here(), "scripts", "diagnostics", "config", "report_config.yaml")
  }
  
  if (!file.exists(config_file)) {
    # Return default config
    return(list(
      report = list(
        n_worst_subjects = 10,
        table_threshold = 50
      )
    ))
  }
  
  yaml::read_yaml(config_file)
}

#' Classify diagnostic value based on thresholds
#' @param value Numeric value to classify
#' @param thresholds List with threshold levels
#' @param metric_type Type of metric ("rhat", "ess_bulk", "ess_tail", "ess_ratio", "divergence_rate", "mcse_ratio", "ebfmi")
#' @return Character string: "PASS", "WARN", or "FAIL"
classify_diagnostic <- function(value, thresholds, metric_type) {
  if (is.na(value) || !is.finite(value)) {
    return("FAIL")
  }
  
  metric_thresholds <- thresholds$thresholds[[metric_type]]
  
  if (metric_type == "rhat") {
    if (value <= metric_thresholds$excellent) return("PASS")
    if (value <= metric_thresholds$acceptable) return("WARN")
    return("FAIL")
  } else if (metric_type == "ess_bulk" || metric_type == "ess_tail") {
    # Use absolute ESS counts (not ratios) for validity assessment
    if (value >= metric_thresholds$good) return("PASS")
    if (value >= metric_thresholds$acceptable) return("WARN")
    return("FAIL")  # < critical threshold
  } else if (metric_type == "ess_ratio") {
    # Ratio is for efficiency info only - don't use for PASS/FAIL
    # Just return informational status
    if (value >= metric_thresholds$efficient) return("PASS")
    if (value >= metric_thresholds$inefficient) return("WARN")
    return("INFO")  # Very inefficient but not invalid
  } else if (metric_type == "divergence_rate") {
    if (value <= metric_thresholds$acceptable) return("PASS")
    if (value <= metric_thresholds$borderline) return("WARN")
    return("FAIL")
  } else if (metric_type == "mcse_ratio") {
    if (value <= metric_thresholds$acceptable) return("PASS")
    if (value <= metric_thresholds$borderline) return("WARN")
    return("FAIL")
  } else if (metric_type == "ebfmi") {
    if (value >= metric_thresholds$acceptable) return("PASS")
    return("FAIL")
  }
  
  return("FAIL")
}

#' Calculate overall status from individual classifications
#' @param classifications Character vector of individual classifications
#' @return Character string: "PASS", "WARN", or "FAIL"
aggregate_status <- function(classifications) {
  if (any(classifications == "FAIL")) {
    return("FAIL")
  }
  if (any(classifications == "WARN")) {
    return("WARN")
  }
  return("PASS")
}

#' Calculate problem score for a subject
#' @param subject_diagnostics List with diagnostic values for one subject
#' @param thresholds List with threshold values
#' @return Numeric problem score
calculate_problem_score <- function(subject_diagnostics, thresholds) {
  score <- 0
  weights <- thresholds$problem_weights
  
  # High Rhat
  if (!is.null(subject_diagnostics$worst_rhat)) {
    if (subject_diagnostics$worst_rhat > thresholds$thresholds$rhat$problematic) {
      score <- score + weights$high_rhat
    }
  }
  
  # Divergences
  if (!is.null(subject_diagnostics$divergence_rate)) {
    if (subject_diagnostics$divergence_rate > thresholds$thresholds$divergence_rate$problematic) {
      score <- score + weights$divergences
    } else if (subject_diagnostics$divergence_rate > thresholds$thresholds$divergence_rate$borderline) {
      score <- score + (weights$divergences * 0.5)
    }
  }
  
  # Low ESS - use absolute counts, not ratios
  if (!is.null(subject_diagnostics$min_ess_bulk)) {
    if (subject_diagnostics$min_ess_bulk < thresholds$thresholds$ess_bulk$critical) {
      score <- score + weights$low_ess
    } else if (subject_diagnostics$min_ess_bulk < thresholds$thresholds$ess_bulk$acceptable) {
      score <- score + (weights$low_ess * 0.5)
    }
  }
  
  if (!is.null(subject_diagnostics$min_ess_tail)) {
    if (subject_diagnostics$min_ess_tail < thresholds$thresholds$ess_tail$critical) {
      score <- score + weights$low_ess
    } else if (subject_diagnostics$min_ess_tail < thresholds$thresholds$ess_tail$acceptable) {
      score <- score + (weights$low_ess * 0.5)
    }
  }
  
  # High MCSE
  if (!is.null(subject_diagnostics$worst_mcse_ratio)) {
    if (subject_diagnostics$worst_mcse_ratio > thresholds$thresholds$mcse_ratio$problematic) {
      score <- score + weights$high_mcse
    } else if (subject_diagnostics$worst_mcse_ratio > thresholds$thresholds$mcse_ratio$borderline) {
      score <- score + (weights$high_mcse * 0.5)
    }
  }
  
  # Low EBFMI
  if (!is.null(subject_diagnostics$ebfmi)) {
    if (subject_diagnostics$ebfmi < thresholds$thresholds$ebfmi$problematic) {
      score <- score + weights$low_ebfmi
    }
  }
  
  return(score)
}

#' Determine fit type from fit object
#' @param fit Fit object to analyze
#' @return Character string: "single", "batch", or "hierarchical"
determine_fit_type <- function(fit) {
  # Check if it's a list of fits (batch)
  if (is.list(fit) && !is.null(names(fit))) {
    # Check if first element is also a fit-like structure
    first_element <- fit[[1]]
    if (is.list(first_element) && any(c("draws", "fit", "diagnostics") %in% names(first_element))) {
      # Check if all elements look like individual fits
      all_fits <- all(sapply(fit, function(x) {
        is.list(x) && any(c("draws", "fit", "diagnostics") %in% names(x))
      }))
      if (all_fits) {
        return("batch")
      }
    }
  }
  
  # Check for hierarchical structure
  # Look for subject-indexed parameters in the draws
  if (!is.null(fit$draws) || !is.null(fit$all_params)) {
    param_names <- if (!is.null(fit$all_params)) fit$all_params else dimnames(fit$draws)[[3]]
    
    # Look for patterns like param[1], param[2], etc. which suggest multiple subjects
    has_indexed_params <- any(grepl("\\[\\d+\\]", param_names))
    
    if (has_indexed_params) {
      # Count unique base parameter names (without indices)
      base_params <- unique(gsub("\\[.*?\\]", "", param_names))
      indexed_params <- param_names[grepl("\\[\\d+\\]", param_names)]
      
      # If we have many indexed versions of the same parameter, it's hierarchical
      if (length(indexed_params) > 10) {  # Arbitrary threshold
        return("hierarchical")
      }
    }
  }
  
  # Default to single
  return("single")
}

#' Extract subject IDs from fit object
#' @param fit Fit object
#' @return Character vector of subject IDs
extract_subject_ids <- function(fit) {
  # First check if subject_list is directly stored
  if (!is.null(fit$subject_list)) {
    return(as.character(fit$subject_list))
  }
  
  # For batch fits, extract from names
  if (is.list(fit) && !is.null(names(fit))) {
    first_element <- fit[[1]]
    if (is.list(first_element) && any(c("draws", "fit") %in% names(first_element))) {
      # It's a batch - names might be subject IDs
      return(names(fit))
    }
  }
  
  # For hierarchical, try to extract from data or metadata
  if (!is.null(fit$data_list) && !is.null(fit$data_list$sid)) {
    return(as.character(fit$data_list$sid))
  }
  
  # If we can't find subject IDs, return NULL
  return(NULL)
}

#' Create diagnostic summary for a single fit
#' @param fit Single fit object
#' @param subject_id Subject ID (optional)
#' @param thresholds Threshold configuration
#' @return List with diagnostic summary
create_single_fit_summary <- function(fit, subject_id = NULL, thresholds = NULL) {
  if (is.null(thresholds)) {
    thresholds <- load_diagnostic_thresholds()
  }
  
  # Use existing diagnostic functions
  source(file.path(here::here(), "scripts", "helpers", "helper_diagnostics.R"))
  
  # Get all parameters or use stored params
  if (!is.null(fit$params)) {
    params <- fit$params
  } else if (!is.null(fit$all_params)) {
    params <- fit$all_params
  } else {
    params <- dimnames(fit$draws)[[3]]
  }
  
  summary <- list(
    subject_id = subject_id,
    n_params = length(params),
    status = "PASS"
  )
  
  # Extract diagnostics
  if (!is.null(fit$diagnostics)) {
    diagnostics <- fit$diagnostics
    
    # Rhat
    if ("rhat" %in% colnames(diagnostics)) {
      rhat_values <- diagnostics[, "rhat"]
      summary$worst_rhat <- max(rhat_values, na.rm = TRUE)
      summary$median_rhat <- median(rhat_values, na.rm = TRUE)
      summary$rhat_status <- classify_diagnostic(summary$worst_rhat, thresholds, "rhat")
    }
    
    # ESS_bulk - use absolute counts for validity, ratio for efficiency info
    if ("ess_bulk" %in% colnames(diagnostics)) {
      ess_bulk <- diagnostics[, "ess_bulk"]
      total_samples <- if (!is.null(fit$tss)) fit$tss else fit$n_iter * fit$n_chains
      ess_bulk_ratio <- ess_bulk / total_samples
      
      # Store both absolute and ratio
      summary$min_ess_bulk <- min(ess_bulk, na.rm = TRUE)
      summary$median_ess_bulk <- median(ess_bulk, na.rm = TRUE)
      summary$min_ess_bulk_ratio <- min(ess_bulk_ratio, na.rm = TRUE)
      summary$median_ess_bulk_ratio <- median(ess_bulk_ratio, na.rm = TRUE)
      
      # Status based on ABSOLUTE ESS, not ratio
      summary$ess_bulk_status <- classify_diagnostic(summary$min_ess_bulk, thresholds, "ess_bulk")
    }
    
    # ESS_tail
    if ("ess_tail" %in% colnames(diagnostics)) {
      ess_tail <- diagnostics[, "ess_tail"]
      total_samples <- if (!is.null(fit$tss)) fit$tss else fit$n_iter * fit$n_chains
      ess_tail_ratio <- ess_tail / total_samples
      
      # Store both absolute and ratio
      summary$min_ess_tail <- min(ess_tail, na.rm = TRUE)
      summary$median_ess_tail <- median(ess_tail, na.rm = TRUE)
      summary$min_ess_tail_ratio <- min(ess_tail_ratio, na.rm = TRUE)
      summary$median_ess_tail_ratio <- median(ess_tail_ratio, na.rm = TRUE)
      
      # Status based on ABSOLUTE ESS, not ratio
      summary$ess_tail_status <- classify_diagnostic(summary$min_ess_tail, thresholds, "ess_tail")
    }
    
    # Overall ESS status is worst of bulk/tail
    if (!is.null(summary$ess_bulk_status) && !is.null(summary$ess_tail_status)) {
      summary$ess_status <- aggregate_status(c(summary$ess_bulk_status, summary$ess_tail_status))
    } else if (!is.null(summary$ess_bulk_status)) {
      summary$ess_status <- summary$ess_bulk_status
    } else if (!is.null(summary$ess_tail_status)) {
      summary$ess_status <- summary$ess_tail_status
    }
  }
  
  # Divergences
  if (!is.null(fit$sampler_diagnostics)) {
    divergences <- sum(fit$sampler_diagnostics[,,"divergent__"], na.rm = TRUE)
    n_iter <- dim(fit$sampler_diagnostics)[1]
    n_chains <- dim(fit$sampler_diagnostics)[2]
    summary$divergence_count <- divergences
    summary$divergence_rate <- divergences / (n_iter * n_chains)
    summary$divergence_status <- classify_diagnostic(summary$divergence_rate, thresholds, "divergence_rate")
  } else if (!is.null(fit$diagnostic_summary$num_divergent)) {
    summary$divergence_count <- sum(fit$diagnostic_summary$num_divergent)
    total_iter <- (fit$n_iter - fit$n_warmup) * fit$n_chains
    summary$divergence_rate <- summary$divergence_count / total_iter
    summary$divergence_status <- classify_diagnostic(summary$divergence_rate, thresholds, "divergence_rate")
  }
  
  # MCSE
  if (!is.null(fit$draws)) {
    draws_matrix <- posterior::as_draws_matrix(fit$draws)
    mcse_values <- apply(draws_matrix, 2, posterior::mcse_mean)
    posterior_sd <- apply(draws_matrix, 2, sd, na.rm = TRUE)
    mcse_ratio <- mcse_values / posterior_sd
    summary$worst_mcse_ratio <- max(mcse_ratio, na.rm = TRUE)
    summary$median_mcse_ratio <- median(mcse_ratio, na.rm = TRUE)
    summary$mcse_status <- classify_diagnostic(summary$worst_mcse_ratio, thresholds, "mcse_ratio")
  }
  
  # EBFMI
  if (!is.null(fit$diagnostic_summary$ebfmi)) {
    summary$ebfmi <- min(fit$diagnostic_summary$ebfmi, na.rm = TRUE)
    summary$ebfmi_status <- classify_diagnostic(summary$ebfmi, thresholds, "ebfmi")
  }
  
  # Calculate overall status
  all_statuses <- c(
    summary$rhat_status,
    summary$ess_status,
    summary$divergence_status,
    summary$mcse_status,
    summary$ebfmi_status
  )
  summary$status <- aggregate_status(all_statuses[!is.na(all_statuses)])
  
  # Calculate problem score
  summary$problem_score <- calculate_problem_score(summary, thresholds)
  
  return(summary)
}

#' Format status with HTML badge
#' @param status Character string: "PASS", "WARN", or "FAIL"
#' @return HTML string with colored badge
format_status_badge <- function(status) {
  colors <- list(
    PASS = "#28a745",  # green
    WARN = "#ffc107",  # yellow
    FAIL = "#dc3545"   # red
  )
  
  color <- colors[[status]]
  if (is.null(color)) color <- "#6c757d"  # gray for unknown
  
  sprintf('<span style="display:inline-block;padding:4px 8px;border-radius:4px;background-color:%s;color:white;font-weight:bold;">%s</span>', 
          color, status)
}

#' Create recommendations based on diagnostic results
#' @param diagnostic_summary Diagnostic summary object
#' @param fit_type Type of fit ("single", "batch", "hierarchical")
#' @return Character vector of recommendations
create_recommendations <- function(diagnostic_summary, fit_type = "single") {
  recommendations <- character()
  
  if (diagnostic_summary$status == "PASS") {
    recommendations <- c(recommendations, 
                        "✓ All diagnostics look good. Proceed with confidence.",
                        "✓ The MCMC sampler has converged and is providing reliable estimates.")
  } else if (diagnostic_summary$status == "WARN") {
    recommendations <- c(recommendations,
                        "⚠ Some diagnostics are borderline. Review carefully before proceeding.",
                        "⚠ Consider running additional iterations or checking model specification.")
  } else {
    recommendations <- c(recommendations,
                        "✗ Critical diagnostic issues detected.",
                        "✗ Do NOT use these results for inference without addressing issues.")
  }
  
  # Specific recommendations based on issues
  if (!is.null(diagnostic_summary$divergence_rate) && 
      diagnostic_summary$divergence_rate > 0.01) {
    recommendations <- c(recommendations,
                        sprintf("• High divergence rate (%.2f%%). Increase adapt_delta to 0.99 or reparameterize model.",
                               diagnostic_summary$divergence_rate * 100))
  }
  
  if (!is.null(diagnostic_summary$worst_rhat) && 
      diagnostic_summary$worst_rhat > 1.1) {
    recommendations <- c(recommendations,
                        sprintf("• High R-hat (%.3f). Run longer chains or investigate multimodality.",
                               diagnostic_summary$worst_rhat))
  }
  
  if (!is.null(diagnostic_summary$min_ess_bulk) && 
      diagnostic_summary$min_ess_bulk < 400) {
    recommendations <- c(recommendations,
                        sprintf("• Low ESS_bulk (%d). Run more iterations for reliable estimates.",
                               as.integer(diagnostic_summary$min_ess_bulk)))
  }
  
  if (!is.null(diagnostic_summary$min_ess_tail) && 
      diagnostic_summary$min_ess_tail < 400) {
    recommendations <- c(recommendations,
                        sprintf("• Low ESS_tail (%d). Run more iterations for reliable quantile estimates.",
                               as.integer(diagnostic_summary$min_ess_tail)))
  }
  
  if (!is.null(diagnostic_summary$ebfmi) && 
      diagnostic_summary$ebfmi < 0.3) {
    recommendations <- c(recommendations,
                        "• Low EBFMI suggests funnel geometry. Consider non-centered parameterization.")
  }
  
  return(recommendations)
}
