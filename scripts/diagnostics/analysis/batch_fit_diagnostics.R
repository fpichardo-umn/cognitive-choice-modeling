# Batch Fit Diagnostics Analysis
# Aggregate diagnostics across multiple subject fits

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
})

# Source required helpers
source(file.path(here::here(), "scripts", "diagnostics", "helpers", "diagnostic_helpers.R"))
source(file.path(here::here(), "scripts", "diagnostics", "analysis", "single_fit_diagnostics.R"))

#' Analyze diagnostics for batch fits
#' @param batch_fit List of fit objects (one per subject)
#' @param thresholds Threshold configuration (optional)
#' @param n_worst Number of worst subjects to highlight (default: 10)
#' @return List with aggregated batch diagnostics
analyze_batch_fits <- function(batch_fit, thresholds = NULL, n_worst = 10) {
  if (is.null(thresholds)) {
    thresholds <- load_diagnostic_thresholds()
  }
  
  # Extract subject IDs
  subject_ids <- names(batch_fit)
  if (is.null(subject_ids)) {
    subject_ids <- as.character(seq_along(batch_fit))
  }
  
  n_subjects <- length(batch_fit)
  
  cat("Analyzing batch of", n_subjects, "subjects...\n")
  
  # Analyze each subject
  subject_analyses <- vector("list", n_subjects)
  subject_summaries <- vector("list", n_subjects)
  
  for (i in seq_along(batch_fit)) {
    if (i %% 50 == 0) {
      cat("  Processed", i, "/", n_subjects, "subjects\n")
    }
    
    tryCatch({
      subject_analyses[[i]] <- analyze_single_fit(batch_fit[[i]], 
                                                   subject_id = subject_ids[i],
                                                   thresholds = thresholds,
                                                   include_plots = FALSE)
      subject_summaries[[i]] <- subject_analyses[[i]]$summary
      subject_summaries[[i]]$subject_id <- subject_ids[i]
    }, error = function(e) {
      warning(sprintf("Failed to analyze subject %s: %s", subject_ids[i], e$message))
      subject_summaries[[i]] <- list(
        subject_id = subject_ids[i],
        status = "ERROR",
        error_message = e$message
      )
    })
  }
  
  # Combine into batch analysis
  batch_analysis <- list(
    n_subjects = n_subjects,
    subject_ids = subject_ids,
    aggregate = aggregate_batch_diagnostics(subject_summaries, thresholds),
    subjects = subject_summaries,
    parameter_summary = aggregate_parameter_diagnostics(subject_analyses),
    problematic_subjects = identify_problematic_subjects_batch(subject_summaries, thresholds, n_worst)
  )
  
  # Overall batch status
  batch_analysis$overall_status <- determine_batch_status(batch_analysis$aggregate)
  batch_analysis$recommendations <- create_batch_recommendations(batch_analysis)
  
  return(batch_analysis)
}

#' Aggregate diagnostics across subjects
#' @param subject_summaries List of subject diagnostic summaries
#' @param thresholds Threshold configuration
#' @return List with aggregated statistics
aggregate_batch_diagnostics <- function(subject_summaries, thresholds) {
  # Extract valid summaries (no errors)
  valid_summaries <- subject_summaries[sapply(subject_summaries, function(x) x$status != "ERROR")]
  n_valid <- length(valid_summaries)
  n_total <- length(subject_summaries)
  
  if (n_valid == 0) {
    return(list(
      n_valid = 0,
      n_total = n_total,
      error = "No valid subjects to analyze"
    ))
  }
  
  # Status counts
  status_counts <- table(sapply(valid_summaries, function(x) x$status))
  
  # R-hat statistics
  rhat_values <- sapply(valid_summaries, function(x) {
    if (!is.null(x$worst_rhat)) x$worst_rhat else NA_real_
  })
  rhat_values <- rhat_values[!is.na(rhat_values)]
  
  # ESS_bulk statistics
  ess_bulk_values <- sapply(valid_summaries, function(x) {
    if (!is.null(x$min_ess_bulk)) x$min_ess_bulk else NA_real_
  })
  ess_bulk_values <- ess_bulk_values[!is.na(ess_bulk_values)]
  
  # ESS_tail statistics
  ess_tail_values <- sapply(valid_summaries, function(x) {
    if (!is.null(x$min_ess_tail)) x$min_ess_tail else NA_real_
  })
  ess_tail_values <- ess_tail_values[!is.na(ess_tail_values)]
  
  # Divergence statistics
  div_rates <- sapply(valid_summaries, function(x) {
    if (!is.null(x$divergence_rate)) x$divergence_rate else NA_real_
  })
  div_rates <- div_rates[!is.na(div_rates)]
  
  # MCSE statistics
  mcse_ratios <- sapply(valid_summaries, function(x) {
    if (!is.null(x$worst_mcse_ratio)) x$worst_mcse_ratio else NA_real_
  })
  mcse_ratios <- mcse_ratios[!is.na(mcse_ratios)]
  
  # EBFMI statistics
  ebfmi_values <- sapply(valid_summaries, function(x) {
    if (!is.null(x$ebfmi)) x$ebfmi else NA_real_
  })
  ebfmi_values <- ebfmi_values[!is.na(ebfmi_values)]
  
  aggregate <- list(
    n_valid = n_valid,
    n_total = n_total,
    n_errors = n_total - n_valid,
    
    status = list(
      pass = as.numeric(status_counts["PASS"]) %||% 0,
      warn = as.numeric(status_counts["WARN"]) %||% 0,
      fail = as.numeric(status_counts["FAIL"]) %||% 0,
      pct_pass = (as.numeric(status_counts["PASS"]) %||% 0) / n_valid * 100,
      pct_warn = (as.numeric(status_counts["WARN"]) %||% 0) / n_valid * 100,
      pct_fail = (as.numeric(status_counts["FAIL"]) %||% 0) / n_valid * 100
    ),
    
    rhat = list(
      min = min(rhat_values, na.rm = TRUE),
      max = max(rhat_values, na.rm = TRUE),
      median = median(rhat_values, na.rm = TRUE),
      mean = mean(rhat_values, na.rm = TRUE),
      q25 = quantile(rhat_values, 0.25, na.rm = TRUE),
      q75 = quantile(rhat_values, 0.75, na.rm = TRUE),
      q90 = quantile(rhat_values, 0.90, na.rm = TRUE),
      n_problematic = sum(rhat_values > thresholds$thresholds$rhat$problematic, na.rm = TRUE),
      pct_problematic = sum(rhat_values > thresholds$thresholds$rhat$problematic, na.rm = TRUE) / length(rhat_values) * 100
    ),
    
    ess_bulk = list(
      min = min(ess_bulk_values, na.rm = TRUE),
      max = max(ess_bulk_values, na.rm = TRUE),
      median = median(ess_bulk_values, na.rm = TRUE),
      mean = mean(ess_bulk_values, na.rm = TRUE),
      q10 = quantile(ess_bulk_values, 0.10, na.rm = TRUE),
      q25 = quantile(ess_bulk_values, 0.25, na.rm = TRUE),
      q75 = quantile(ess_bulk_values, 0.75, na.rm = TRUE),
      n_critical = sum(ess_bulk_values < thresholds$thresholds$ess_bulk$critical, na.rm = TRUE),
      pct_critical = sum(ess_bulk_values < thresholds$thresholds$ess_bulk$critical, na.rm = TRUE) / length(ess_bulk_values) * 100
    ),
    
    ess_tail = list(
      min = min(ess_tail_values, na.rm = TRUE),
      max = max(ess_tail_values, na.rm = TRUE),
      median = median(ess_tail_values, na.rm = TRUE),
      mean = mean(ess_tail_values, na.rm = TRUE),
      q10 = quantile(ess_tail_values, 0.10, na.rm = TRUE),
      q25 = quantile(ess_tail_values, 0.25, na.rm = TRUE),
      q75 = quantile(ess_tail_values, 0.75, na.rm = TRUE),
      n_critical = sum(ess_tail_values < thresholds$thresholds$ess_tail$critical, na.rm = TRUE),
      pct_critical = sum(ess_tail_values < thresholds$thresholds$ess_tail$critical, na.rm = TRUE) / length(ess_tail_values) * 100
    ),
    
    divergences = list(
      min_rate = min(div_rates, na.rm = TRUE),
      max_rate = max(div_rates, na.rm = TRUE),
      median_rate = median(div_rates, na.rm = TRUE),
      mean_rate = mean(div_rates, na.rm = TRUE),
      n_with_divergences = sum(div_rates > 0, na.rm = TRUE),
      pct_with_divergences = sum(div_rates > 0, na.rm = TRUE) / length(div_rates) * 100,
      n_problematic = sum(div_rates > thresholds$thresholds$divergence_rate$problematic, na.rm = TRUE),
      pct_problematic = sum(div_rates > thresholds$thresholds$divergence_rate$problematic, na.rm = TRUE) / length(div_rates) * 100
    ),
    
    mcse = list(
      min_ratio = min(mcse_ratios, na.rm = TRUE),
      max_ratio = max(mcse_ratios, na.rm = TRUE),
      median_ratio = median(mcse_ratios, na.rm = TRUE),
      mean_ratio = mean(mcse_ratios, na.rm = TRUE),
      q90 = quantile(mcse_ratios, 0.90, na.rm = TRUE),
      n_problematic = sum(mcse_ratios > thresholds$thresholds$mcse_ratio$problematic, na.rm = TRUE),
      pct_problematic = sum(mcse_ratios > thresholds$thresholds$mcse_ratio$problematic, na.rm = TRUE) / length(mcse_ratios) * 100
    )
  )
  
  if (length(ebfmi_values) > 0) {
    aggregate$ebfmi <- list(
      min = min(ebfmi_values, na.rm = TRUE),
      median = median(ebfmi_values, na.rm = TRUE),
      mean = mean(ebfmi_values, na.rm = TRUE),
      n_problematic = sum(ebfmi_values < thresholds$thresholds$ebfmi$problematic, na.rm = TRUE),
      pct_problematic = sum(ebfmi_values < thresholds$thresholds$ebfmi$problematic, na.rm = TRUE) / length(ebfmi_values) * 100
    )
  }
  
  return(aggregate)
}

#' Aggregate parameter diagnostics across subjects
#' @param subject_analyses List of subject analysis results
#' @return List with parameter-level aggregation
aggregate_parameter_diagnostics <- function(subject_analyses) {
  # Get valid analyses
  valid_analyses <- subject_analyses[!sapply(subject_analyses, is.null)]
  
  if (length(valid_analyses) == 0) {
    return(NULL)
  }
  
  # Get all parameter names (assuming consistent across subjects)
  all_params <- valid_analyses[[1]]$parameters$parameter
  n_params <- length(all_params)
  n_subjects <- length(valid_analyses)
  
  # Initialize aggregation matrices
  rhat_matrix <- matrix(NA_real_, nrow = n_subjects, ncol = n_params)
  ess_bulk_matrix <- matrix(NA_real_, nrow = n_subjects, ncol = n_params)
  ess_tail_matrix <- matrix(NA_real_, nrow = n_subjects, ncol = n_params)
  mcse_matrix <- matrix(NA_real_, nrow = n_subjects, ncol = n_params)
  
  colnames(rhat_matrix) <- all_params
  colnames(ess_bulk_matrix) <- all_params
  colnames(ess_tail_matrix) <- all_params
  colnames(mcse_matrix) <- all_params
  
  # Fill matrices
  for (i in seq_along(valid_analyses)) {
    param_df <- valid_analyses[[i]]$parameters
    rhat_matrix[i, ] <- param_df$rhat
    ess_bulk_matrix[i, ] <- param_df$ess_bulk  # Use absolute counts
    ess_tail_matrix[i, ] <- param_df$ess_tail  # Use absolute counts
    mcse_matrix[i, ] <- param_df$mcse_ratio
  }
  
  # Calculate per-parameter statistics
  param_summary <- data.frame(
    parameter = all_params,
    rhat_mean = colMeans(rhat_matrix, na.rm = TRUE),
    rhat_median = apply(rhat_matrix, 2, median, na.rm = TRUE),
    rhat_max = apply(rhat_matrix, 2, max, na.rm = TRUE),
    rhat_n_problematic = colSums(rhat_matrix > 1.1, na.rm = TRUE),
    ess_bulk_mean = colMeans(ess_bulk_matrix, na.rm = TRUE),
    ess_bulk_median = apply(ess_bulk_matrix, 2, median, na.rm = TRUE),
    ess_bulk_min = apply(ess_bulk_matrix, 2, min, na.rm = TRUE),
    ess_bulk_n_critical = colSums(ess_bulk_matrix < 100, na.rm = TRUE),
    ess_tail_mean = colMeans(ess_tail_matrix, na.rm = TRUE),
    ess_tail_median = apply(ess_tail_matrix, 2, median, na.rm = TRUE),
    ess_tail_min = apply(ess_tail_matrix, 2, min, na.rm = TRUE),
    ess_tail_n_critical = colSums(ess_tail_matrix < 100, na.rm = TRUE),
    mcse_mean = colMeans(mcse_matrix, na.rm = TRUE),
    mcse_median = apply(mcse_matrix, 2, median, na.rm = TRUE),
    mcse_max = apply(mcse_matrix, 2, max, na.rm = TRUE),
    mcse_n_problematic = colSums(mcse_matrix > 0.15, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Identify consistently problematic parameters
  param_summary$is_problematic <- (
    param_summary$rhat_n_problematic > n_subjects * 0.1 |
    param_summary$ess_bulk_n_critical > n_subjects * 0.1 |
    param_summary$ess_tail_n_critical > n_subjects * 0.1 |
    param_summary$mcse_n_problematic > n_subjects * 0.1
  )
  
  return(list(
    summary = param_summary,
    problematic_params = param_summary$parameter[param_summary$is_problematic],
    n_problematic_params = sum(param_summary$is_problematic)
  ))
}

#' Identify problematic subjects in batch
#' @param subject_summaries List of subject summaries
#' @param thresholds Threshold configuration
#' @param n_worst Number of worst subjects to return
#' @return Data frame with problematic subjects
identify_problematic_subjects_batch <- function(subject_summaries, thresholds, n_worst) {
  # Filter out errors
  valid_summaries <- subject_summaries[sapply(subject_summaries, function(x) x$status != "ERROR")]
  
  if (length(valid_summaries) == 0) {
    return(NULL)
  }
  
  # Create subject table
  subject_table <- data.frame(
    subject_id = sapply(valid_summaries, function(x) x$subject_id),
    status = sapply(valid_summaries, function(x) x$status),
    worst_rhat = sapply(valid_summaries, function(x) x$worst_rhat %||% NA_real_),
    median_rhat = sapply(valid_summaries, function(x) x$median_rhat %||% NA_real_),
    min_ess_bulk = sapply(valid_summaries, function(x) x$min_ess_bulk %||% NA_real_),
    min_ess_tail = sapply(valid_summaries, function(x) x$min_ess_tail %||% NA_real_),
    min_ess_bulk_ratio = sapply(valid_summaries, function(x) x$min_ess_bulk_ratio %||% NA_real_),
    min_ess_tail_ratio = sapply(valid_summaries, function(x) x$min_ess_tail_ratio %||% NA_real_),
    divergence_count = sapply(valid_summaries, function(x) x$divergence_count %||% 0),
    divergence_rate = sapply(valid_summaries, function(x) x$divergence_rate %||% 0),
    worst_mcse_ratio = sapply(valid_summaries, function(x) x$worst_mcse_ratio %||% NA_real_),
    median_mcse_ratio = sapply(valid_summaries, function(x) x$median_mcse_ratio %||% NA_real_),
    ebfmi = sapply(valid_summaries, function(x) x$ebfmi %||% NA_real_),
    problem_score = sapply(valid_summaries, function(x) x$problem_score %||% 0),
    stringsAsFactors = FALSE
  )
  
  # Sort by problem score
  subject_table <- subject_table[order(subject_table$problem_score, decreasing = TRUE), ]
  
  # Identify problematic parameters for each subject (if available)
  # This would require more detailed parameter info which we'll skip for now
  
  # Return top N worst
  n_to_return <- min(n_worst, nrow(subject_table))
  worst_subjects <- subject_table[1:n_to_return, ]
  
  # Add all subjects table for potential CSV export
  return(list(
    worst = worst_subjects,
    all = subject_table,
    n_problematic = sum(subject_table$status %in% c("WARN", "FAIL")),
    pct_problematic = sum(subject_table$status %in% c("WARN", "FAIL")) / nrow(subject_table) * 100
  ))
}

#' Determine overall batch status
#' @param aggregate Aggregated batch diagnostics
#' @return Character string: "PASS", "WARN", or "FAIL"
determine_batch_status <- function(aggregate) {
  # Handle missing percentages
  pct_fail <- aggregate$status$pct_fail %||% 0
  pct_warn <- aggregate$status$pct_warn %||% 0
  
  # If more than 20% failed, batch fails
  if (!is.na(pct_fail) && pct_fail > 20) {
    return("FAIL")
  }
  
  # If more than 50% have warnings or failures, batch fails
  if (!is.na(pct_fail) && !is.na(pct_warn) && (pct_warn + pct_fail) > 50) {
    return("FAIL")
  }
  
  # If more than 10% failed or more than 30% warned, batch warns
  if ((!is.na(pct_fail) && pct_fail > 10) || (!is.na(pct_warn) && pct_warn > 30)) {
    return("WARN")
  }
  
  # Otherwise pass
  return("PASS")
}

#' Create recommendations for batch analysis
#' @param batch_analysis Batch analysis results
#' @return Character vector of recommendations
create_batch_recommendations <- function(batch_analysis) {
  recommendations <- character()
  agg <- batch_analysis$aggregate
  
  # Overall status
  if (batch_analysis$overall_status == "PASS") {
    recommendations <- c(recommendations,
                        "✓ Batch diagnostics look good overall.",
                        sprintf("✓ %d/%d subjects (%.1f%%) passed all diagnostics.",
                               agg$status$pass, agg$n_valid, agg$status$pct_pass))
  } else if (batch_analysis$overall_status == "WARN") {
    recommendations <- c(recommendations,
                        "⚠ Some subjects have diagnostic issues.",
                        "⚠ Review problematic subjects before proceeding with analysis.")
  } else {
    recommendations <- c(recommendations,
                        "✗ Significant diagnostic issues across batch.",
                        "✗ Substantial portion of subjects have convergence problems.")
  }
  
  # Specific issues
  if (!is.null(agg$status$pct_fail) && !is.na(agg$status$pct_fail) && agg$status$pct_fail > 5) {
    recommendations <- c(recommendations,
                        sprintf("• %.1f%% of subjects failed diagnostics. Consider refitting these subjects.",
                               agg$status$pct_fail))
  }
  
  if (!is.null(agg$rhat$pct_problematic) && !is.na(agg$rhat$pct_problematic) && agg$rhat$pct_problematic > 10) {
    recommendations <- c(recommendations,
                        sprintf("• %.1f%% of subjects have high R-hat values (>1.1).",
                               agg$rhat$pct_problematic))
  }
  
  if (!is.null(agg$divergences$pct_problematic) && !is.na(agg$divergences$pct_problematic) && agg$divergences$pct_problematic > 5) {
    recommendations <- c(recommendations,
                        sprintf("• %.1f%% of subjects have high divergence rates.",
                               agg$divergences$pct_problematic),
                        "  Consider increasing adapt_delta or reparameterizing model.")
  }
  
  if (!is.null(agg$ess_bulk) && !is.null(agg$ess_bulk$pct_critical) && 
      !is.na(agg$ess_bulk$pct_critical) && agg$ess_bulk$pct_critical > 10) {
    recommendations <- c(recommendations,
                        sprintf("• %.1f%% of subjects have critically low ESS_bulk (<100).",
                               agg$ess_bulk$pct_critical),
                        "  Run more iterations for reliable estimates.")
  }
  
  if (!is.null(agg$ess_tail) && !is.null(agg$ess_tail$pct_critical) && 
      !is.na(agg$ess_tail$pct_critical) && agg$ess_tail$pct_critical > 10) {
    recommendations <- c(recommendations,
                        sprintf("• %.1f%% of subjects have critically low ESS_tail (<100).",
                               agg$ess_tail$pct_critical),
                        "  Run more iterations for reliable quantile estimates.")
  }
  
  # Parameter-specific recommendations
  if (!is.null(batch_analysis$parameter_summary) && 
      batch_analysis$parameter_summary$n_problematic_params > 0) {
    recommendations <- c(recommendations,
                        sprintf("• %d parameters consistently problematic across subjects:",
                               batch_analysis$parameter_summary$n_problematic_params),
                        sprintf("  %s", paste(head(batch_analysis$parameter_summary$problematic_params, 5), 
                                             collapse = ", ")))
  }
  
  return(recommendations)
}

#' Create subject diagnostic table for export
#' @param batch_analysis Batch analysis results
#' @param include_all Whether to include all subjects or just problematic (default: TRUE)
#' @return Data frame ready for CSV export
create_subject_table_for_export <- function(batch_analysis, include_all = TRUE) {
  if (include_all) {
    return(batch_analysis$problematic_subjects$all)
  } else {
    prob_subjects <- batch_analysis$problematic_subjects$all
    prob_subjects <- prob_subjects[prob_subjects$status %in% c("WARN", "FAIL"), ]
    return(prob_subjects)
  }
}
