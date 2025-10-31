# Problematic Subjects Identification and Ranking
# Cross-cutting analysis for identifying subjects with diagnostic issues

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})

# Source required helpers
source(file.path(here::here(), "scripts", "diagnostics", "helpers", "diagnostic_helpers.R"))

#' Identify and rank problematic subjects
#' @param diagnostic_results Results from batch or hierarchical analysis
#' @param thresholds Threshold configuration (optional)
#' @param n_worst Number of worst subjects to return (default: 10)
#' @param fit_type Type of fit analysis ("batch" or "hierarchical")
#' @return List with problematic subject analysis
identify_problematic_subjects <- function(diagnostic_results, thresholds = NULL, 
                                         n_worst = 10, fit_type = "batch") {
  if (is.null(thresholds)) {
    thresholds <- load_diagnostic_thresholds()
  }
  
  # Extract subject-level information based on fit type
  if (fit_type == "batch") {
    subject_info <- extract_subjects_from_batch(diagnostic_results)
  } else if (fit_type == "hierarchical") {
    subject_info <- extract_subjects_from_hierarchical(diagnostic_results)
  } else {
    stop("fit_type must be 'batch' or 'hierarchical'")
  }
  
  if (is.null(subject_info) || nrow(subject_info) == 0) {
    return(NULL)
  }
  
  # Calculate comprehensive problem scores
  subject_info <- calculate_comprehensive_scores(subject_info, thresholds)
  
  # Rank subjects by problem severity
  subject_info <- subject_info[order(subject_info$problem_score, decreasing = TRUE), ]
  
  # Categorize problems
  problem_categories <- categorize_problems(subject_info, thresholds)
  
  # Get top N worst subjects
  n_return <- min(n_worst, nrow(subject_info))
  worst_subjects <- subject_info[1:n_return, ]
  
  # Identify subjects with specific problem types
  specific_problems <- identify_specific_problem_types(subject_info, thresholds)
  
  # Create summary
  summary <- list(
    n_total_subjects = nrow(subject_info),
    n_problematic = sum(subject_info$status %in% c("WARN", "FAIL")),
    pct_problematic = sum(subject_info$status %in% c("WARN", "FAIL")) / nrow(subject_info) * 100,
    n_fail = sum(subject_info$status == "FAIL"),
    n_warn = sum(subject_info$status == "WARN"),
    n_pass = sum(subject_info$status == "PASS"),
    problem_categories = problem_categories,
    specific_problems = specific_problems
  )
  
  return(list(
    summary = summary,
    all_subjects = subject_info,
    worst_subjects = worst_subjects,
    subjects_by_problem = specific_problems
  ))
}

#' Extract subject information from batch analysis
#' @param batch_results Batch analysis results
#' @return Data frame with subject information
extract_subjects_from_batch <- function(batch_results) {
  if (is.null(batch_results$problematic_subjects$all)) {
    return(NULL)
  }
  
  return(batch_results$problematic_subjects$all)
}

#' Extract subject information from hierarchical analysis
#' @param hier_results Hierarchical analysis results
#' @return Data frame with subject information
extract_subjects_from_hierarchical <- function(hier_results) {
  if (is.null(hier_results$subject_analysis$subject_summaries)) {
    return(NULL)
  }
  
  # Convert list of summaries to data frame
  summaries <- hier_results$subject_analysis$subject_summaries
  
  subject_df <- data.frame(
    subject_id = sapply(summaries, function(x) x$subject_id %||% NA_character_),
    status = sapply(summaries, function(x) x$status %||% NA_character_),
    worst_rhat = sapply(summaries, function(x) x$worst_rhat %||% NA_real_),
    median_rhat = sapply(summaries, function(x) x$median_rhat %||% NA_real_),
    min_ess_ratio = sapply(summaries, function(x) x$min_ess_ratio %||% NA_real_),
    problem_score = sapply(summaries, function(x) x$problem_score %||% 0),
    stringsAsFactors = FALSE
  )
  
  return(subject_df)
}

#' Calculate comprehensive problem scores
#' @param subject_info Data frame with subject information
#' @param thresholds Threshold configuration
#' @return Data frame with updated problem scores
calculate_comprehensive_scores <- function(subject_info, thresholds) {
  # If problem_score already calculated, we're done
  if ("problem_score" %in% colnames(subject_info) && all(!is.na(subject_info$problem_score))) {
    return(subject_info)
  }
  
  # Otherwise calculate
  weights <- thresholds$problem_weights
  
  subject_info$problem_score <- 0
  
  # R-hat contribution
  if ("worst_rhat" %in% colnames(subject_info)) {
    subject_info$problem_score <- subject_info$problem_score +
      ifelse(subject_info$worst_rhat > thresholds$thresholds$rhat$problematic, weights$high_rhat, 0)
  }
  
  # Divergence contribution
  if ("divergence_rate" %in% colnames(subject_info)) {
    subject_info$problem_score <- subject_info$problem_score +
      ifelse(subject_info$divergence_rate > thresholds$thresholds$divergence_rate$problematic, 
             weights$divergences,
             ifelse(subject_info$divergence_rate > thresholds$thresholds$divergence_rate$borderline,
                   weights$divergences * 0.5, 0))
  }
  
  # ESS contribution
  if ("min_ess_ratio" %in% colnames(subject_info)) {
    subject_info$problem_score <- subject_info$problem_score +
      ifelse(subject_info$min_ess_ratio < thresholds$thresholds$ess_ratio$problematic,
             weights$low_ess,
             ifelse(subject_info$min_ess_ratio < thresholds$thresholds$ess_ratio$acceptable,
                   weights$low_ess * 0.5, 0))
  }
  
  # MCSE contribution
  if ("worst_mcse_ratio" %in% colnames(subject_info)) {
    subject_info$problem_score <- subject_info$problem_score +
      ifelse(subject_info$worst_mcse_ratio > thresholds$thresholds$mcse_ratio$problematic,
             weights$high_mcse,
             ifelse(subject_info$worst_mcse_ratio > thresholds$thresholds$mcse_ratio$borderline,
                   weights$high_mcse * 0.5, 0))
  }
  
  # EBFMI contribution
  if ("ebfmi" %in% colnames(subject_info)) {
    subject_info$problem_score <- subject_info$problem_score +
      ifelse(!is.na(subject_info$ebfmi) & subject_info$ebfmi < thresholds$thresholds$ebfmi$problematic,
             weights$low_ebfmi, 0)
  }
  
  return(subject_info)
}

#' Categorize problems across subjects
#' @param subject_info Data frame with subject information
#' @param thresholds Threshold configuration
#' @return List with problem category counts
categorize_problems <- function(subject_info, thresholds) {
  categories <- list(
    high_rhat = 0,
    low_ess = 0,
    high_divergences = 0,
    high_mcse = 0,
    low_ebfmi = 0,
    multiple_issues = 0
  )
  
  for (i in 1:nrow(subject_info)) {
    row <- subject_info[i, ]
    issue_count <- 0
    
    if ("worst_rhat" %in% colnames(subject_info) && 
        !is.na(row$worst_rhat) && 
        row$worst_rhat > thresholds$thresholds$rhat$problematic) {
      categories$high_rhat <- categories$high_rhat + 1
      issue_count <- issue_count + 1
    }
    
    if ("min_ess_ratio" %in% colnames(subject_info) && 
        !is.na(row$min_ess_ratio) && 
        row$min_ess_ratio < thresholds$thresholds$ess_ratio$problematic) {
      categories$low_ess <- categories$low_ess + 1
      issue_count <- issue_count + 1
    }
    
    if ("divergence_rate" %in% colnames(subject_info) && 
        !is.na(row$divergence_rate) && 
        row$divergence_rate > thresholds$thresholds$divergence_rate$problematic) {
      categories$high_divergences <- categories$high_divergences + 1
      issue_count <- issue_count + 1
    }
    
    if ("worst_mcse_ratio" %in% colnames(subject_info) && 
        !is.na(row$worst_mcse_ratio) && 
        row$worst_mcse_ratio > thresholds$thresholds$mcse_ratio$problematic) {
      categories$high_mcse <- categories$high_mcse + 1
      issue_count <- issue_count + 1
    }
    
    if ("ebfmi" %in% colnames(subject_info) && 
        !is.na(row$ebfmi) && 
        row$ebfmi < thresholds$thresholds$ebfmi$problematic) {
      categories$low_ebfmi <- categories$low_ebfmi + 1
      issue_count <- issue_count + 1
    }
    
    if (issue_count > 1) {
      categories$multiple_issues <- categories$multiple_issues + 1
    }
  }
  
  return(categories)
}

#' Identify subjects with specific problem types
#' @param subject_info Data frame with subject information
#' @param thresholds Threshold configuration
#' @return List of subject IDs by problem type
identify_specific_problem_types <- function(subject_info, thresholds) {
  problems <- list(
    high_rhat_only = character(),
    low_ess_only = character(),
    high_divergences_only = character(),
    multiple_issues = character(),
    convergence_issues = character(),
    sampling_issues = character()
  )
  
  for (i in 1:nrow(subject_info)) {
    row <- subject_info[i, ]
    subj_id <- row$subject_id
    
    # Check individual issues
    has_rhat <- "worst_rhat" %in% colnames(subject_info) && 
      !is.na(row$worst_rhat) && 
      row$worst_rhat > thresholds$thresholds$rhat$problematic
    
    has_ess <- "min_ess_ratio" %in% colnames(subject_info) && 
      !is.na(row$min_ess_ratio) && 
      row$min_ess_ratio < thresholds$thresholds$ess_ratio$problematic
    
    has_div <- "divergence_rate" %in% colnames(subject_info) && 
      !is.na(row$divergence_rate) && 
      row$divergence_rate > thresholds$thresholds$divergence_rate$problematic
    
    issue_count <- sum(c(has_rhat, has_ess, has_div))
    
    # Categorize
    if (issue_count == 0) {
      next
    } else if (issue_count == 1) {
      if (has_rhat) problems$high_rhat_only <- c(problems$high_rhat_only, subj_id)
      if (has_ess) problems$low_ess_only <- c(problems$low_ess_only, subj_id)
      if (has_div) problems$high_divergences_only <- c(problems$high_divergences_only, subj_id)
    } else {
      problems$multiple_issues <- c(problems$multiple_issues, subj_id)
    }
    
    # Broader categories
    if (has_rhat || has_ess) {
      problems$convergence_issues <- c(problems$convergence_issues, subj_id)
    }
    
    if (has_div) {
      problems$sampling_issues <- c(problems$sampling_issues, subj_id)
    }
  }
  
  return(problems)
}

#' Create comparison table for cross-model analysis
#' @param problematic_subjects_list List of problematic subjects from different models
#' @param model_names Names of models
#' @return Data frame with cross-model comparison
create_cross_model_comparison <- function(problematic_subjects_list, model_names) {
  if (length(problematic_subjects_list) != length(model_names)) {
    stop("Number of models must match number of model names")
  }
  
  # Get all unique subject IDs
  all_subjects <- unique(unlist(lapply(problematic_subjects_list, function(x) x$all_subjects$subject_id)))
  
  # Create comparison data frame
  comparison <- data.frame(
    subject_id = all_subjects,
    stringsAsFactors = FALSE
  )
  
  # Add columns for each model
  for (i in seq_along(model_names)) {
    model_name <- model_names[i]
    model_subjects <- problematic_subjects_list[[i]]$all_subjects
    
    # Match subjects
    match_idx <- match(comparison$subject_id, model_subjects$subject_id)
    
    comparison[[paste0(model_name, "_status")]] <- model_subjects$status[match_idx]
    comparison[[paste0(model_name, "_score")]] <- model_subjects$problem_score[match_idx]
    comparison[[paste0(model_name, "_rhat")]] <- if ("worst_rhat" %in% colnames(model_subjects)) {
      model_subjects$worst_rhat[match_idx]
    } else NA_real_
  }
  
  # Calculate cross-model statistics
  score_cols <- grep("_score$", colnames(comparison), value = TRUE)
  comparison$mean_score <- rowMeans(comparison[, score_cols], na.rm = TRUE)
  comparison$n_models_problematic <- rowSums(!is.na(comparison[, score_cols]) & 
                                              comparison[, score_cols] > 0)
  
  # Sort by mean problem score
  comparison <- comparison[order(comparison$mean_score, decreasing = TRUE), ]
  
  return(comparison)
}

#' Generate problem summary text
#' @param problematic_analysis Results from identify_problematic_subjects
#' @return Character vector with summary text
generate_problem_summary_text <- function(problematic_analysis) {
  summary <- problematic_analysis$summary
  
  text <- c(
    sprintf("Total subjects analyzed: %d", summary$n_total_subjects),
    sprintf("Subjects with issues: %d (%.1f%%)", summary$n_problematic, summary$pct_problematic),
    sprintf("  - Failed: %d", summary$n_fail),
    sprintf("  - Warned: %d", summary$n_warn),
    sprintf("  - Passed: %d", summary$n_pass),
    "",
    "Problem breakdown:",
    sprintf("  - High R-hat: %d subjects", summary$problem_categories$high_rhat),
    sprintf("  - Low ESS: %d subjects", summary$problem_categories$low_ess),
    sprintf("  - High divergences: %d subjects", summary$problem_categories$high_divergences),
    sprintf("  - Multiple issues: %d subjects", summary$problem_categories$multiple_issues)
  )
  
  return(text)
}

#' Export problematic subjects to CSV
#' @param problematic_analysis Results from identify_problematic_subjects
#' @param output_file Path to output CSV file
#' @param include_all Whether to include all subjects or just problematic (default: TRUE)
#' @return TRUE if export successful
export_problematic_subjects_csv <- function(problematic_analysis, output_file, include_all = TRUE) {
  if (include_all) {
    data_to_export <- problematic_analysis$all_subjects
  } else {
    data_to_export <- problematic_analysis$all_subjects[
      problematic_analysis$all_subjects$status %in% c("WARN", "FAIL"), 
    ]
  }
  
  tryCatch({
    write.csv(data_to_export, output_file, row.names = FALSE)
    return(TRUE)
  }, error = function(e) {
    warning("Failed to export CSV: ", e$message)
    return(FALSE)
  })
}
