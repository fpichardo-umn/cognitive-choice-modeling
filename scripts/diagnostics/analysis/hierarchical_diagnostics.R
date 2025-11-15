# Hierarchical Fit Diagnostics Analysis
# Analyze diagnostics for hierarchical models (fit-based or real hierarchical)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(posterior)
})

# Source required helpers
source(file.path(here::here(), "scripts", "diagnostics", "helpers", "diagnostic_helpers.R"))
source(file.path(here::here(), "scripts", "diagnostics", "analysis", "single_fit_diagnostics.R"))

#' Analyze hierarchical fit diagnostics
#' @param hier_fit Hierarchical fit object
#' @param thresholds Threshold configuration (optional)
#' @param real_hier Whether to treat as real hierarchical (default: FALSE)
#' @param n_worst Number of worst subjects to highlight (default: 10)
#' @return List with hierarchical diagnostics
analyze_hierarchical_fit <- function(hier_fit, thresholds = NULL, real_hier = FALSE, n_worst = 10) {
  if (is.null(thresholds)) {
    thresholds <- load_diagnostic_thresholds()
  }
  
  # Identify subject-level and group-level parameters
  param_structure <- identify_hierarchical_parameters(hier_fit)
  
  cat("Analyzing hierarchical fit with", param_structure$n_subjects, "subjects...\n")
  
  # Always analyze subject-level parameters
  subject_analysis <- analyze_subject_level_parameters(hier_fit, param_structure, thresholds, n_worst)
  
  hier_analysis <- list(
    fit_type = if (real_hier) "real_hierarchical" else "fit_based_hierarchical",
    n_subjects = param_structure$n_subjects,
    subject_analysis = subject_analysis,
    basic_checks = check_basic_hierarchical_health(hier_fit, param_structure, thresholds)
  )
  
  # Always analyze group-level parameters (lp__, mu_pr, sigma, mu_*)
  hier_analysis$group_analysis <- analyze_group_level_parameters(hier_fit, param_structure, thresholds)
  
  # If real hierarchical, do full shrinkage analysis
  if (real_hier) {
    hier_analysis$shrinkage_analysis <- analyze_shrinkage_full(hier_fit, param_structure)
  } else {
    # For fit-based, just basic shrinkage check
    hier_analysis$shrinkage_check <- check_shrinkage_basic(hier_fit, param_structure)
  }
  
  # Overall assessment
  hier_analysis$overall_status <- determine_hierarchical_status(hier_analysis)
  hier_analysis$recommendations <- create_hierarchical_recommendations(hier_analysis, real_hier)
  
  return(hier_analysis)
}

#' Identify hierarchical parameter structure
#' @param hier_fit Hierarchical fit object
#' @return List with parameter structure info
identify_hierarchical_parameters <- function(hier_fit) {
  all_params <- if (!is.null(hier_fit$all_params)) {
    hier_fit$all_params
  } else {
    dimnames(hier_fit$draws)[[3]]
  }
  
  # Keep lp__ for diagnostics
  has_lp <- "lp__" %in% all_params
  
  # Identify indexed parameters (e.g., alpha[1], mu_pr[2])
  indexed_params <- all_params[grepl("\\[\\d+\\]", all_params)]
  
  # CRITICAL: Separate hyperparameters (indexed by param) from subject params (indexed by subject)
  # mu_pr[*] and sigma[*] are group-level hyperparameters, NOT subject parameters
  hyperparam_indexed <- indexed_params[grepl("^(mu_pr|sigma)\\[", indexed_params)]
  subject_params <- setdiff(indexed_params, hyperparam_indexed)
  
  # Extract unique base parameter names for subject-level only
  base_params <- unique(gsub("\\[\\d+\\]", "", subject_params))
  
  # Infer number of subjects from subject parameter indices
  if (length(subject_params) > 0) {
    all_indices <- as.numeric(gsub(".*\\[(\\d+)\\].*", "\\1", subject_params))
    n_subjects <- max(all_indices, na.rm = TRUE)
  } else {
    n_subjects <- 1
  }
  
  # Identify generated mu quantities (e.g., mu_gain, mu_loss)
  # These are scalar (non-indexed) and follow pattern mu_<param_name>
  non_indexed_params <- all_params[!grepl("\\[\\d+\\]", all_params)]
  non_indexed_params <- non_indexed_params[non_indexed_params != "lp__"]
  generated_mus <- non_indexed_params[grepl("^mu_[a-z_]+$", non_indexed_params)]
  
  # Group-level parameters include:
  # 1. Indexed hyperparameters: mu_pr[*], sigma[*]
  # 2. Generated quantities: mu_gain, mu_loss, etc.
  # 3. lp__ (if present)
  group_params <- c(hyperparam_indexed, generated_mus)
  if (has_lp) {
    group_params <- c("lp__", group_params)
  }
  
  structure <- list(
    n_subjects = n_subjects,
    base_params = base_params,
    n_base_params = length(base_params),
    subject_params = subject_params,
    group_params = group_params,
    hyperparam_indexed = hyperparam_indexed,
    generated_mus = generated_mus,
    has_lp = has_lp,
    n_group_params = length(group_params),
    all_params = all_params
  )
  
  return(structure)
}

#' Analyze subject-level parameters in hierarchical model
#' @param hier_fit Hierarchical fit object
#' @param param_structure Parameter structure from identify_hierarchical_parameters
#' @param thresholds Threshold configuration
#' @param n_worst Number of worst subjects
#' @return List with subject-level analysis
analyze_subject_level_parameters <- function(hier_fit, param_structure, thresholds, n_worst) {
  subject_params <- param_structure$subject_params
  n_subjects <- param_structure$n_subjects
  
  # Extract diagnostics for subject parameters
  if (!is.null(hier_fit$diagnostics)) {
    # Get all params from fit (keep lp__ in the list)
    fit_all_params <- if (!is.null(hier_fit$all_params)) {
      hier_fit$all_params
    } else {
      rownames(hier_fit$diagnostics)
    }
    
    # Get indices for subject parameters
    subj_param_indices <- which(fit_all_params %in% subject_params)
    
    if (length(subj_param_indices) > 0) {
      subj_diagnostics <- hier_fit$diagnostics[subj_param_indices, , drop = FALSE]
      
      # Organize by subject
      subject_summaries <- vector("list", n_subjects)
      
      for (i in 1:n_subjects) {
        # Get parameters for this subject
        subj_i_params <- subject_params[grepl(sprintf("\\[%d\\]", i), subject_params)]
        subj_i_indices <- which(fit_all_params %in% subj_i_params)
        
        if (length(subj_i_indices) > 0) {
          subj_i_diag <- subj_diagnostics[rownames(subj_diagnostics) %in% subj_i_params, , drop = FALSE]
          
          # Calculate summary for this subject
          subject_summaries[[i]] <- list(
            subject_id = if (!is.null(hier_fit$subject_list)) hier_fit$subject_list[i] else as.character(i),
            worst_rhat = max(subj_i_diag[, "rhat"], na.rm = TRUE),
            median_rhat = median(subj_i_diag[, "rhat"], na.rm = TRUE),
            # Store both absolute ESS and ratio
            min_ess_bulk = if ("ess_bulk" %in% colnames(subj_i_diag)) {
              min(subj_i_diag[, "ess_bulk"], na.rm = TRUE)
            } else NA_real_,
            min_ess_bulk_ratio = if ("ess_bulk" %in% colnames(subj_i_diag)) {
              total_samples <- hier_fit$n_iter * hier_fit$n_chains
              min(subj_i_diag[, "ess_bulk"] / total_samples, na.rm = TRUE)
            } else NA_real_,
            min_ess_tail = if ("ess_tail" %in% colnames(subj_i_diag)) {
              min(subj_i_diag[, "ess_tail"], na.rm = TRUE)
            } else NA_real_,
            status = "PASS"  # Will be updated based on classifications
          )
          
          # Classify status based on ABSOLUTE ESS, not ratio
          rhat_status <- classify_diagnostic(subject_summaries[[i]]$worst_rhat, thresholds, "rhat")
          ess_bulk_status <- if (!is.na(subject_summaries[[i]]$min_ess_bulk)) {
            classify_diagnostic(subject_summaries[[i]]$min_ess_bulk, thresholds, "ess_bulk")
          } else "PASS"
          ess_tail_status <- if (!is.na(subject_summaries[[i]]$min_ess_tail)) {
            classify_diagnostic(subject_summaries[[i]]$min_ess_tail, thresholds, "ess_tail")
          } else "PASS"
          ess_status <- aggregate_status(c(ess_bulk_status, ess_tail_status))
          
          subject_summaries[[i]]$status <- aggregate_status(c(rhat_status, ess_status))
          subject_summaries[[i]]$problem_score <- calculate_problem_score(subject_summaries[[i]], thresholds)
        }
      }
      
      # Aggregate across subjects
      aggregate_subj <- list(
        n_subjects = n_subjects,
        n_params_per_subject = length(param_structure$base_params),
        
        rhat = list(
          min = min(sapply(subject_summaries, function(x) x$worst_rhat), na.rm = TRUE),
          max = max(sapply(subject_summaries, function(x) x$worst_rhat), na.rm = TRUE),
          median = median(sapply(subject_summaries, function(x) x$worst_rhat), na.rm = TRUE),
          n_problematic = sum(sapply(subject_summaries, function(x) x$worst_rhat > thresholds$thresholds$rhat$problematic), na.rm = TRUE)
        ),
        
        ess_bulk = list(
          min = min(sapply(subject_summaries, function(x) x$min_ess_bulk), na.rm = TRUE),
          median = median(sapply(subject_summaries, function(x) x$min_ess_bulk), na.rm = TRUE),
          n_critical = sum(sapply(subject_summaries, function(x) x$min_ess_bulk < thresholds$thresholds$ess_bulk$critical), na.rm = TRUE)
        ),
        
        ess_tail = list(
          min = min(sapply(subject_summaries, function(x) x$min_ess_tail), na.rm = TRUE),
          median = median(sapply(subject_summaries, function(x) x$min_ess_tail), na.rm = TRUE),
          n_critical = sum(sapply(subject_summaries, function(x) x$min_ess_tail < thresholds$thresholds$ess_tail$critical), na.rm = TRUE)
        ),
        
        status = list(
          pass = sum(sapply(subject_summaries, function(x) x$status == "PASS")),
          warn = sum(sapply(subject_summaries, function(x) x$status == "WARN")),
          fail = sum(sapply(subject_summaries, function(x) x$status == "FAIL"))
        )
      )
      
      # Identify worst subjects
      subject_scores <- sapply(subject_summaries, function(x) x$problem_score)
      worst_indices <- order(subject_scores, decreasing = TRUE)[1:min(n_worst, n_subjects)]
      
      return(list(
        aggregate = aggregate_subj,
        subject_summaries = subject_summaries,
        worst_subjects = subject_summaries[worst_indices]
      ))
    }
  }
  
  return(NULL)
}

#' Basic hierarchical health checks
#' @param hier_fit Hierarchical fit object
#' @param param_structure Parameter structure
#' @param thresholds Threshold configuration
#' @return List with basic health checks
check_basic_hierarchical_health <- function(hier_fit, param_structure, thresholds) {
  checks <- list()
  
  # Check if hyperparameters are not degenerate
  if (length(param_structure$group_params) > 0 && !is.null(hier_fit$diagnostics)) {
    # Get all params from fit (keep lp__ now)
    fit_all_params <- if (!is.null(hier_fit$all_params)) {
      hier_fit$all_params
    } else {
      rownames(hier_fit$diagnostics)
    }
    
    group_param_indices <- which(fit_all_params %in% param_structure$group_params)
    
    if (length(group_param_indices) > 0) {
      group_diagnostics <- hier_fit$diagnostics[group_param_indices, , drop = FALSE]
      
      checks$hyperparameters <- list(
        converged = all(group_diagnostics[, "rhat"] < thresholds$thresholds$rhat$acceptable, na.rm = TRUE),
        worst_rhat = max(group_diagnostics[, "rhat"], na.rm = TRUE),
        status = if (all(group_diagnostics[, "rhat"] < thresholds$thresholds$rhat$acceptable, na.rm = TRUE)) "PASS" else "FAIL"
      )
    }
  }
  
  # Check overall divergences
  if (!is.null(hier_fit$diagnostic_summary$num_divergent)) {
    total_div <- sum(hier_fit$diagnostic_summary$num_divergent)
    # n_iter is already post-warmup samples per chain
    total_trans <- hier_fit$n_iter * hier_fit$n_chains
    div_rate <- total_div / total_trans
    
    checks$divergences <- list(
      count = total_div,
      rate = div_rate,
      status = classify_diagnostic(div_rate, thresholds, "divergence_rate")
    )
  }
  
  return(checks)
}

#' Basic shrinkage check for fit-based hierarchical
#' @param hier_fit Hierarchical fit object
#' @param param_structure Parameter structure
#' @return List with basic shrinkage assessment
check_shrinkage_basic <- function(hier_fit, param_structure) {
  # This is a simplified check - just verify that we're getting reasonable variation
  # across subjects (not all shrunk to same value, not completely unregularized)
  
  if (length(param_structure$subject_params) == 0) {
    return(NULL)
  }
  
  # For each base parameter, check variance across subjects
  shrinkage_info <- list()
  
  for (base_param in param_structure$base_params) {
    # Get all subject versions of this parameter
    param_pattern <- sprintf("%s\\[\\d+\\]", base_param)
    subj_params <- param_structure$subject_params[grepl(param_pattern, param_structure$subject_params)]
    
    if (length(subj_params) > 1 && !is.null(hier_fit$draws)) {
      # Get posterior means for each subject
      param_indices <- which(dimnames(hier_fit$draws)[[3]] %in% subj_params)
      param_draws <- hier_fit$draws[, , param_indices, drop = FALSE]
      param_means <- apply(param_draws, 3, mean, na.rm = TRUE)
      
      shrinkage_info[[base_param]] <- list(
        variance = var(param_means, na.rm = TRUE),
        range = diff(range(param_means, na.rm = TRUE)),
        cv = sd(param_means, na.rm = TRUE) / abs(mean(param_means, na.rm = TRUE))
      )
    }
  }
  
  return(list(
    working = length(shrinkage_info) > 0,
    parameter_info = shrinkage_info
  ))
}

#' Analyze group-level parameters (for real hierarchical)
#' @param hier_fit Hierarchical fit object
#' @param param_structure Parameter structure
#' @param thresholds Threshold configuration
#' @return List with group-level analysis
analyze_group_level_parameters <- function(hier_fit, param_structure, thresholds) {
  if (length(param_structure$group_params) == 0) {
    return(NULL)
  }
  
  # Get all params from fit (keep lp__ now)
  fit_all_params <- if (!is.null(hier_fit$all_params)) {
    hier_fit$all_params
  } else if (!is.null(hier_fit$diagnostics)) {
    rownames(hier_fit$diagnostics)
  } else {
    NULL
  }
  
  group_param_indices <- which(fit_all_params %in% param_structure$group_params)
  
  if (length(group_param_indices) > 0 && !is.null(hier_fit$diagnostics)) {
    group_diagnostics <- hier_fit$diagnostics[group_param_indices, , drop = FALSE]
    
    # Categorize parameters for clearer reporting
    param_categories <- character(length(param_structure$group_params))
    names(param_categories) <- param_structure$group_params
    
    for (p in param_structure$group_params) {
      if (p == "lp__") {
        param_categories[p] <- "log_posterior"
      } else if (grepl("^mu_pr\\[", p)) {
        param_categories[p] <- "hyper_mean"
      } else if (grepl("^sigma\\[", p)) {
        param_categories[p] <- "hyper_sd"
      } else if (grepl("^mu_", p)) {
        param_categories[p] <- "generated_mu"
      } else {
        param_categories[p] <- "other"
      }
    }
    
    # Extract R-hat, ESS, etc for group parameters
    group_analysis <- list(
      n_params = length(param_structure$group_params),
      parameters = param_structure$group_params,
      param_categories = param_categories,
      
      rhat = list(
        values = group_diagnostics[, "rhat"],
        max = max(group_diagnostics[, "rhat"], na.rm = TRUE),
        all_converged = all(group_diagnostics[, "rhat"] < thresholds$thresholds$rhat$acceptable, na.rm = TRUE)
      ),
      
      ess_bulk = if ("ess_bulk" %in% colnames(group_diagnostics)) {
        list(
          values = group_diagnostics[, "ess_bulk"],
          min = min(group_diagnostics[, "ess_bulk"], na.rm = TRUE)
        )
      } else NULL,
      
      ess_tail = if ("ess_tail" %in% colnames(group_diagnostics)) {
        list(
          values = group_diagnostics[, "ess_tail"],
          min = min(group_diagnostics[, "ess_tail"], na.rm = TRUE)
        )
      } else NULL
    )
    
    # Overall group-level status based on absolute ESS thresholds
    rhat_ok <- group_analysis$rhat$all_converged
    ess_bulk_ok <- is.null(group_analysis$ess_bulk) || 
      group_analysis$ess_bulk$min >= thresholds$thresholds$ess_bulk$acceptable
    ess_tail_ok <- is.null(group_analysis$ess_tail) || 
      group_analysis$ess_tail$min >= thresholds$thresholds$ess_tail$acceptable
    
    group_analysis$status <- if (rhat_ok && ess_bulk_ok && ess_tail_ok) {
      "PASS"
    } else if (group_analysis$rhat$max < thresholds$thresholds$rhat$problematic) {
      "WARN"
    } else {
      "FAIL"
    }
    
    return(group_analysis)
  }
  
  return(NULL)
}

#' Full shrinkage analysis for real hierarchical
#' @param hier_fit Hierarchical fit object
#' @param param_structure Parameter structure
#' @return List with detailed shrinkage analysis
analyze_shrinkage_full <- function(hier_fit, param_structure) {
  # This would include:
  # - Shrinkage factor calculations
  # - Comparison of individual vs group estimates
  # - Between-subject variation assessment
  # - Effective number of independent subjects
  
  # For now, simplified version
  shrinkage <- check_shrinkage_basic(hier_fit, param_structure)
  
  # Could add: comparison to non-hierarchical fit if available
  # Could add: formal shrinkage factor calculations
  
  return(shrinkage)
}

#' Determine overall hierarchical model status
#' @param hier_analysis Hierarchical analysis results
#' @return Character string: "PASS", "WARN", or "FAIL"
determine_hierarchical_status <- function(hier_analysis) {
  # Check basic health
  if (!is.null(hier_analysis$basic_checks$hyperparameters) &&
      hier_analysis$basic_checks$hyperparameters$status == "FAIL") {
    return("FAIL")
  }
  
  if (!is.null(hier_analysis$basic_checks$divergences) &&
      hier_analysis$basic_checks$divergences$status == "FAIL") {
    return("FAIL")
  }
  
  # Check subject-level
  if (!is.null(hier_analysis$subject_analysis$aggregate$status)) {
    fail_pct <- hier_analysis$subject_analysis$aggregate$status$fail / 
      hier_analysis$n_subjects * 100
    warn_pct <- hier_analysis$subject_analysis$aggregate$status$warn / 
      hier_analysis$n_subjects * 100
    
    if (fail_pct > 20) return("FAIL")
    if (fail_pct > 10 || warn_pct > 30) return("WARN")
  }
  
  # Check group-level if real hierarchical
  if (!is.null(hier_analysis$group_analysis) &&
      hier_analysis$group_analysis$status %in% c("FAIL", "WARN")) {
    return(hier_analysis$group_analysis$status)
  }
  
  return("PASS")
}

#' Create recommendations for hierarchical model
#' @param hier_analysis Hierarchical analysis results
#' @param real_hier Whether this is real hierarchical
#' @return Character vector of recommendations
create_hierarchical_recommendations <- function(hier_analysis, real_hier) {
  recommendations <- character()
  
  # Overall status
  if (hier_analysis$overall_status == "PASS") {
    recommendations <- c(recommendations,
                        "✓ Hierarchical model diagnostics look good.",
                        "✓ Subject-level parameters have converged.")
  } else if (hier_analysis$overall_status == "WARN") {
    recommendations <- c(recommendations,
                        "⚠ Some diagnostic issues in hierarchical model.",
                        "⚠ Review specific subjects or parameters before proceeding.")
  } else {
    recommendations <- c(recommendations,
                        "✗ Significant diagnostic issues in hierarchical model.",
                        "✗ Model may need reparameterization or longer sampling.")
  }
  
  # Hyperparameter checks
  if (!is.null(hier_analysis$basic_checks$hyperparameters) &&
      hier_analysis$basic_checks$hyperparameters$status != "PASS") {
    recommendations <- c(recommendations,
                        sprintf("• Hyperparameters not fully converged (worst R-hat: %.3f).",
                               hier_analysis$basic_checks$hyperparameters$worst_rhat))
  }
  
  # Subject-level issues
  if (!is.null(hier_analysis$subject_analysis$aggregate$status)) {
    fail_count <- hier_analysis$subject_analysis$aggregate$status$fail
    if (fail_count > 0) {
      recommendations <- c(recommendations,
                          sprintf("• %d subjects have failed diagnostics.", fail_count))
    }
  }
  
  # Group-level issues (if real hierarchical)
  if (real_hier && !is.null(hier_analysis$group_analysis) &&
      hier_analysis$group_analysis$status != "PASS") {
    recommendations <- c(recommendations,
                        "• Group-level parameters have convergence issues.",
                        "  This affects population-level inference.")
  }
  
  # Shrinkage assessment
  if (!is.null(hier_analysis$shrinkage_check) && !hier_analysis$shrinkage_check$working) {
    recommendations <- c(recommendations,
                        "⚠ Warning: Hierarchical structure may not be working as intended.",
                        "  Check that there is meaningful variation across subjects.")
  }
  
  return(recommendations)
}
