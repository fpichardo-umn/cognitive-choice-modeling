# Helper functions for model fitting
# These functions handle model fitting, checkpointing, and related utilities

# Load required libraries
suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(here)
})

# Import required helper functions
source(file.path(here::here(), "scripts", "helpers", "helper_dirs.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_functions.R"))

# ============================================================================
# Constants for adaptive iteration feature
# ============================================================================

# Hardcoded problematic thresholds (trigger restart)
PROBLEMATIC_THRESHOLDS <- list(
  min_ess = 100,           # Fatal ESS < 100
  max_divergence_pct = 1.0,  # Fatal 1% divergence rate
  max_rhat = 1.10,           # Fatal R-hat (hardcoded)
  inefficient_ess_threshold = 200 # Rhat > target AND ESS < this triggers restart
)

# Hardcoded restart settings
RESTART_WARMUP_INCREMENT <- 500
RESTART_ADAPT_DELTA_INCREMENT <- 0.025
RESTART_MAX_ADAPT_DELTA <- 0.99 # Kept original per user request
RESTART_MAX_TREEDEPTH_INCREMENT <- 2
RESTART_CHAIN_INCREMENT <- 2
MAX_RESTARTS <- 2

# Default acceptable thresholds (user can override)
# Kept original per user request
DEFAULT_DIAG_THRESHOLDS <- list(
  rhat = 1.01,
  ess_bulk = 400,
  ess_tail = 400
)

# ============================================================================
# Helper functions for adaptive iteration feature
# ============================================================================

#' Check diagnostics for all parameters
#' @param draws Posterior draws object
#' @param diag_thresholds List with rhat, ess_bulk, ess_tail thresholds
#' @return List with all_pass flag and worst values
check_diagnostics <- function(draws, diag_thresholds) {
  # Calculate summary statistics
  summary_stats <- posterior::summarize_draws(draws)
  
  # Extract diagnostics for all parameters
  rhats <- summary_stats$rhat
  ess_bulks <- summary_stats$ess_bulk
  ess_tails <- summary_stats$ess_tail
  
  # Remove NA values (can happen with some parameters)
  rhats <- rhats[!is.na(rhats)]
  ess_bulks <- ess_bulks[!is.na(ess_bulks)]
  ess_tails <- ess_tails[!is.na(ess_tails)]
  
  # Check if all pass
  rhat_pass <- all(rhats < diag_thresholds$rhat)
  ess_bulk_pass <- all(ess_bulks >= diag_thresholds$ess_bulk)
  ess_tail_pass <- all(ess_tails >= diag_thresholds$ess_tail)
  
  all_pass <- rhat_pass && ess_bulk_pass && ess_tail_pass
  
  return(list(
    all_pass = all_pass,
    worst_rhat = if (length(rhats) > 0) max(rhats) else 0,
    worst_ess_bulk = if (length(ess_bulks) > 0) min(ess_bulks) else Inf,
    worst_ess_tail = if (length(ess_tails) > 0) min(ess_tails) else Inf
  ))
}

#' Detect problem type based on diagnostics
#' @param diag_results Results from check_diagnostics()
#' @param num_divergent Total number of divergent transitions
#' @param total_samples Total number of post-warmup samples
#' @param diag_thresholds The *target* thresholds (used for Rhat check)
#' @return String problem type: "pass", "divergence", "rhat_fatal", "ess_low", 
#'         "rhat_and_low_ess", "extension_needed"
get_problem_type <- function(diag_results, num_divergent, total_samples, diag_thresholds) {
  
  # 1. Check Fatal Divergence
  if (total_samples > 0) {
    divergence_pct <- (num_divergent / total_samples) * 100
    if (divergence_pct > PROBLEMATIC_THRESHOLDS$max_divergence_pct) {
      return("divergence")
    }
  } else if (num_divergent > 0) {
    return("divergence") # Divergences with 0 samples is bad
  }
  
  # 2. Check Fatal R-hat
  if (diag_results$worst_rhat >= PROBLEMATIC_THRESHOLDS$max_rhat) {
    return("rhat_fatal")
  }
  
  # 3. Check Fatal ESS
  if (diag_results$worst_ess_bulk < PROBLEMATIC_THRESHOLDS$min_ess ||
      diag_results$worst_ess_tail < PROBLEMATIC_THRESHOLDS$min_ess) {
    return("ess_low")
  }
  
  # 4. Check "Inefficient Convergence" (Rhat > target AND ESS < 200)
  is_rhat_fail <- diag_results$worst_rhat >= diag_thresholds$rhat
  is_ess_inefficient <- diag_results$worst_ess_bulk < PROBLEMATIC_THRESHOLDS$inefficient_ess_threshold ||
    diag_results$worst_ess_tail < PROBLEMATIC_THRESHOLDS$inefficient_ess_threshold
  
  if (is_rhat_fail && is_ess_inefficient) {
    return("rhat_and_low_ess")
  }
  
  # 5. Check for any other failure (e.g., Rhat > 1.01 but ESS > 200, or ESS < 1000)
  if (!diag_results$all_pass) {
    return("extension_needed")
  }
  
  # 6. Success
  return("pass")
}


#' Calculate escalated settings for a restart attempt
#' @param original_warmup Original warmup iterations
#' @param original_adapt_delta Original adapt_delta value
#' @param original_max_treedepth Original max_treedepth value
#' @param original_n_chains Original number of chains
#' @param restart_count Current restart attempt number (1 or 2)
#' @param problem_type The string code for the problem detected
#' @return List with new warmup, adapt_delta, max_treedepth, n_chains
calculate_restart_settings <- function(original_warmup, original_adapt_delta, 
                                       original_max_treedepth, original_n_chains, 
                                       restart_count, problem_type) {
  
  # Start with original settings
  new_warmup <- original_warmup
  new_adapt_delta <- original_adapt_delta
  new_max_treedepth <- original_max_treedepth
  new_n_chains <- original_n_chains
  
  # Apply escalations based on the problem
  if (problem_type == "divergence") {
    new_warmup <- original_warmup + (restart_count * RESTART_WARMUP_INCREMENT)
    new_adapt_delta <- min(original_adapt_delta + (restart_count * RESTART_ADAPT_DELTA_INCREMENT), 
                           RESTART_MAX_ADAPT_DELTA)
    new_max_treedepth <- original_max_treedepth + (restart_count * RESTART_MAX_TREEDEPTH_INCREMENT)
    
  } else if (problem_type == "rhat_fatal") {
    new_warmup <- original_warmup + (restart_count * RESTART_WARMUP_INCREMENT)
    
  } else if (problem_type == "ess_low") {
    new_warmup <- original_warmup + (restart_count * RESTART_WARMUP_INCREMENT)
    new_n_chains <- original_n_chains + (restart_count * RESTART_CHAIN_INCREMENT)
    
  } else if (problem_type == "rhat_and_low_ess") {
    # Fix is same as "ess_low": more warmup and more chains
    new_warmup <- original_warmup + (restart_count * RESTART_WARMUP_INCREMENT)
    new_n_chains <- original_n_chains + (restart_count * RESTART_CHAIN_INCREMENT)
  }
  
  return(list(
    warmup = new_warmup,
    adapt_delta = new_adapt_delta,
    max_treedepth = new_max_treedepth,
    n_chains = new_n_chains
  ))
}

#' Remove duplicate cohort directories from path
#' @param file_path Character string with file path
#' @param cohort Character string with cohort name
#' @return Character string with corrected path
fix_duplicate_cohort_path <- function(file_path, cohort) {
  # Create the duplicate pattern to look for
  duplicate_pattern <- paste0("/", cohort, "/", cohort, "/")
  single_pattern <- paste0("/", cohort, "/")
  
  # Replace the duplicate with single occurrence
  if (grepl(duplicate_pattern, file_path)) {
    fixed_path <- gsub(duplicate_pattern, single_pattern, file_path)
    cat("Fixed duplicate cohort path:", file_path, "->", fixed_path, "\n")
    return(fixed_path)
  }
  
  return(file_path)
}

#' Fit and save a model
#' @param task Character string specifying the task name
#' @param group_type Character string specifying the group type (e.g., "sing", "hier")
#' @param model_name Character string specifying the model name
#' @param model_type Character string specifying the model type (e.g., "fit", "postpc")
#' @param data_list List of data to fit the model to
#' @param n_subs Integer number of subjects
#' @param n_trials Integer number of trials
#' @param n_warmup Integer number of warmup iterations (MCMC only)
#' @param n_iter Integer number of sampling iterations (MCMC only)
#' @param n_chains Integer number of chains (MCMC) or paths (Pathfinder)
#' @param adapt_delta Numeric adaptation delta parameter (MCMC only)
#' @param max_treedepth Integer maximum tree depth (MCMC only)
#' @param model_params Character vector of model parameters
#' @param dry_run Boolean to perform a dry run without actual sampling
#' @param checkpoint_interval Integer number of iterations between checkpoints (MCMC only)
#' @param output_dir Optional character string with output directory path
#' @param emp_bayes Boolean indicating whether this is an empirical Bayes model
#' @param informative_priors Optional list of informative priors
#' @param init_params Optional list of parameter initializations
#' @param fitting_method Character string specifying method: "mcmc" or "pathfinder" (default: "mcmc")
#' @param pf_num_paths Integer number of pathfinder paths (default: 4)
#' @param pf_draws Integer number of final draws after PSIS (default: 1000)
#' @param pf_single_path_draws Integer draws per single path (default: 250)
#' @param pf_psis_resample Boolean whether to use PSIS resampling (default: TRUE)
#' @param pf_calculate_lp Boolean whether to calculate log probability (default: TRUE)
#' @param min_iter Integer minimum post-warmup iterations for adaptive fitting (default: NULL)
#' @param max_iter Integer maximum post-warmup iterations for adaptive fitting (default: NULL)
#' @param iter_increment Integer iterations to add at each diagnostic check (default: 1000)
#' @param diag_thresholds List with rhat, ess_bulk, ess_tail thresholds (default: NULL uses defaults)
#' @param enable_adaptive_iter Boolean to enable adaptive iteration feature (default: TRUE)
#' @return List with fitted model results or NULL for dry run
fit_and_save_model <- function(task, cohort, ses, group_type, model_name, model_type, data_list, 
                               n_subs, n_trials, n_warmup, n_iter, n_chains, adapt_delta, max_treedepth, 
                               model_params, dry_run = FALSE, checkpoint_interval = 1000, 
                               output_dir = NULL, emp_bayes = FALSE, informative_priors = NULL,
                               init_params = NULL, cohort_sub_dir = TRUE, model_status = NULL,
                               is_simulation = FALSE, index = NULL, data_filt_list = NULL,
                               fitting_method = "mcmc", pf_num_paths = 4, pf_draws = 1000,
                               pf_single_path_draws = 250, pf_psis_resample = TRUE, pf_calculate_lp = TRUE,
                               min_iter = NULL, max_iter = NULL, iter_increment = 1000,
                               diag_thresholds = NULL, enable_adaptive_iter = TRUE) {
  
  # Create the model string and get the model path
  model_str <- paste(task, group_type, model_name, sep="_")
  model_path <- get_model_file_path(task, group_type, model_name, model_type, status = model_status)
  
  # Load the Stan model
  tryCatch({
    stanmodel_arg <- cmdstan_model(exe_file = model_path)
  }, error = function(e) {
    stop(paste("Failed to load model from:", model_path, "\nError:", e$message))
  })
  
  # Basic validation for init_params
  if (!is.null(init_params)) {
    if (!is.list(init_params)) {
      stop("init_params must be a list or list of lists")
    }
    
    # Check if it's a list of lists
    if (is.list(init_params[[1]])) {
      if (length(init_params) != n_chains) {
        warning("Number of initialization lists (", length(init_params), 
                ") does not match number of chains (", n_chains, "). Adjusting init_params.")
        # Replicate the first init list to match n_chains
        init_params <- replicate(n_chains, init_params[[1]], simplify = FALSE)
      }
    } else {
      # If single list, replicate it for each chain
      init_params <- replicate(n_chains, init_params, simplify = FALSE)
    }
  }
  
  # Validation and setup for adaptive iteration feature
  use_adaptive_iter <- FALSE
  if (enable_adaptive_iter && !is.null(min_iter) && !is.null(max_iter)) {
    # Validate parameters
    if (min_iter > max_iter) {
      stop("min_iter (", min_iter, ") must be <= max_iter (", max_iter, ")")
    }
    if (min_iter < 1) {
      stop("min_iter must be >= 1")
    }
    if (iter_increment < 1) {
      stop("iter_increment must be >= 1")
    }
    
    # Set up diagnostic thresholds (Use passed-in or default)
    if (is.null(diag_thresholds)) {
      diag_thresholds <- DEFAULT_DIAG_THRESHOLDS
      cat("Using default diagnostic thresholds.\n")
    } else {
      # Validate provided thresholds
      if (!all(c("rhat", "ess_bulk", "ess_tail") %in% names(diag_thresholds))) {
        stop("diag_thresholds must contain rhat, ess_bulk, and ess_tail")
      }
      cat("Using user-supplied diagnostic thresholds.\n")
    }
    
    use_adaptive_iter <- TRUE
    cat("Adaptive iteration enabled: min=", min_iter, ", max=", max_iter, 
        ", increment=", iter_increment, "\n", sep="")
    cat(sprintf("  Target thresholds: Rhat < %.2f, ESS_Bulk >= %d, ESS_Tail >= %d\n",
                diag_thresholds$rhat, diag_thresholds$ess_bulk, diag_thresholds$ess_tail))
  } else if (enable_adaptive_iter && (is.null(min_iter) || is.null(max_iter))) {
    # Adaptive requested but missing parameters - use old behavior
    cat("Adaptive iteration requested but min_iter/max_iter not provided. Using fixed iterations.\n")
  }
  
  # Determine output directory and file paths
  if (is.null(output_dir)) {
    if (emp_bayes) {
      # Empirical Bayes individual fits
      output_dir <- get_empbayes_output_dir(task, "individual", cohort)
    } else if (is_simulation) {
      # Fits to simulated data (for parameter recovery)
      output_dir <- get_validation_output_dir(task, "parameter_recovery", "fits")
      # Could add simulation run ID subfolder here
      output_dir <- file.path(output_dir, paste0("sim-", format(Sys.Date(), "%Y%m%d")))
    } else {
      # Regular fits organized by cohort/session
      output_dir <- get_fits_output_dir(task, cohort, ses)
    }
  }
  
  # Ensure output directory exists
  ensure_dir_exists(output_dir)
  
  # Generate output and checkpoint file paths
  output_file <- get_output_file_path(task, cohort, group_type, model_name, model_type, 
                                      emp_bayes, if (length(data_list$sid) == 1) data_list$sid else NULL,
                                      index, ses = ses, 
                                      output_dir = output_dir, cohort_sub_dir)
  output_file <- fix_duplicate_cohort_path(output_file, cohort)
  checkpoint_file <- paste0(tools::file_path_sans_ext(output_file), "_checkpoint.rds")
  
  # Dry run only prints information without fitting
  if (dry_run) {
    cat("Dry run for model:", model_str, "\n")
    cat("Fitting method:", fitting_method, "\n")
    cat("Data list:\n")
    print(str(data_list))
    cat("Parameters:\n")
    cat("n_subs =", n_subs, "\n")
    if (fitting_method == "mcmc") {
      cat("n_warmup =", n_warmup, "\n")
      cat("n_iter =", n_iter, "\n")
      cat("n_chains =", n_chains, "\n")
      cat("adapt_delta =", adapt_delta, "\n")
      cat("max_treedepth =", max_treedepth, "\n")
    } else if (fitting_method == "pathfinder") {
      cat("num_paths =", pf_num_paths, "\n")
      cat("draws =", pf_draws, "\n")
      cat("single_path_draws =", pf_single_path_draws, "\n")
    }
    if (!is.null(informative_priors)) {
      cat("Using informative priors\n")
    }
    if (!is.null(init_params)) {
      cat("Using custom parameter initialization\n")
    }
    cat("\n", "File to be saved as:", output_file, "\n")
    cat("\n", "Output directory exists:", dir.exists(output_dir), "\n")
    return(NULL)
  }
  
  cat("Fitting model:", model_str, "\n")
  cat("Fitting method:", fitting_method, "\n")
  
  # Add informative priors to data_list if provided
  if (!is.null(informative_priors)) {
    data_list$pr_mu <- unname(sapply(model_params, function(param) informative_priors[[param]][1]))
    data_list$pr_sigma <- unname(sapply(model_params, function(param) informative_priors[[param]][2]))
  }
  
  # Branch based on fitting method
  if (fitting_method == "pathfinder") {
    # === Pathfinder Fitting (No Checkpointing) ===
    
    fit_result <- run_pathfinder_and_process(
      stanmodel_arg, data_list, pf_num_paths, pf_draws, pf_single_path_draws,
      pf_psis_resample, pf_calculate_lp, model_str, task, n_subs, model_params, 
      output_file, init_params
    )
    
  } else if (fitting_method != "mcmc") {
    stop("Invalid fitting_method: must be 'mcmc' or 'pathfinder'")
  } else if (use_adaptive_iter) {
    # === MCMC with Adaptive Iteration ===
    
    restart_count <- 0
    restart_attempts <- list()
    converged <- FALSE
    
    # Original settings (for restart escalation)
    original_warmup <- n_warmup
    original_adapt_delta <- adapt_delta
    original_max_treedepth <- max_treedepth
    original_n_chains <- n_chains
    # This is the "double post warmup" chunk from your ESS rule
    initial_post_warmup_chunk <- min_iter 
    
    # Current settings for this attempt
    current_warmup <- original_warmup
    current_adapt_delta <- original_adapt_delta
    current_max_treedepth <- original_max_treedepth
    current_n_chains <- original_n_chains
    
    # Outer restart loop
    while (restart_count <= MAX_RESTARTS && !converged) {
      
      # Reset state for this restart attempt
      current_iter <- 0
      accumulated_samples <- list()
      accumulated_diagnostics <- list()
      step_size <- NULL
      inv_metric <- NULL
      warmup_done <- FALSE
      extension_count <- 0
      
      # Calculate total iterations needed (warmup + post-warmup)
      # min_iter and max_iter represent POST-WARMUP samples
      target_total_min <- current_warmup + min_iter
      target_total_max <- current_warmup + max_iter
      
      # Run to min_iter first
      cat(sprintf("Running warmup (%d) + minimum post-warmup iterations (%d) = %d total...\n", 
                  current_warmup, min_iter, target_total_min))
      target_iter <- target_total_min
      
      # Checkpoint loop to reach target_iter
      while (current_iter < target_iter) {
        remaining_iter <- min(checkpoint_interval, target_iter - current_iter)
        
        if (!warmup_done) {
          sampling_result <- run_initial_sampling(
            stanmodel_arg, data_list, current_warmup, remaining_iter, current_n_chains,
            current_adapt_delta, current_max_treedepth, init_params
          )
          
          fit <- sampling_result$fit
          new_samples <- sampling_result$new_samples
          current_iter <- sampling_result$current_iter
          
          step_size <- fit$metadata()$step_size_adaptation
          inv_metric <- fit$inv_metric()
          warmup_done <- TRUE
        } else {
          sampling_result <- continue_sampling(
            stanmodel_arg, data_list, remaining_iter, current_n_chains,
            current_adapt_delta, current_max_treedepth,
            accumulated_samples, step_size, inv_metric
          )
          
          fit <- sampling_result$fit
          new_samples <- sampling_result$new_samples
          current_iter <- current_iter + sampling_result$iterations_added
        }
        
        accumulated_samples <- c(accumulated_samples, list(new_samples))
        new_diagnostics <- fit$sampler_diagnostics()
        accumulated_diagnostics <- c(accumulated_diagnostics, list(new_diagnostics))
        
        save_checkpoint(checkpoint_file, current_iter, accumulated_samples,
                        accumulated_diagnostics, step_size, inv_metric, warmup_done)
        
        cat("Checkpoint saved at iteration", current_iter, "\n")
      }
      
      # Reached min_iter - check diagnostics
      cat(sprintf("Iteration %d complete. Checking diagnostics...\n", current_iter))
      
      all_samples <- do.call(posterior::bind_draws, c(accumulated_samples, along = "iteration"))
      all_diagnostics <- do.call(posterior::bind_draws, c(accumulated_diagnostics, along = "iteration"))
      
      diag_results <- check_diagnostics(all_samples, diag_thresholds)
      num_divergent <- colSums(all_diagnostics[, , "divergent__"])
      
      # Calculate total POST-WARMUP samples for divergence rate
      total_post_warmup_samples <- current_n_chains * (current_iter - current_warmup)
      
      # Get the problem type
      problem_type <- get_problem_type(diag_results, sum(num_divergent), 
                                       total_post_warmup_samples, diag_thresholds)
      
      cat(sprintf("  Worst Rhat: %.3f, Min ESS: %.0f, Problem: %s\n",
                  diag_results$worst_rhat,
                  min(c(diag_results$worst_ess_bulk, diag_results$worst_ess_tail), na.rm = TRUE),
                  problem_type))
      
      # --- New Decision Tree ---
      
      # 1. Did we pass?
      if (problem_type == "pass") {
        cat("  All diagnostics passed.\n")
        converged <- TRUE
        break  # Exit restart loop
      }
      
      # 2. Is it a restart-level problem?
      if (problem_type %in% c("divergence", "rhat_fatal", "ess_low", "rhat_and_low_ess")) {
        cat(sprintf("  Fundamental problem '%s' detected.\n", problem_type))
        
        if (restart_count < MAX_RESTARTS) {
          # Save this failed attempt
          restart_attempts[[restart_count + 1]] <- list(
            warmup = current_warmup,
            adapt_delta = current_adapt_delta,
            max_treedepth = current_max_treedepth,
            n_chains = current_n_chains,
            iterations = current_iter,
            worst_rhat = diag_results$worst_rhat,
            min_ess = min(diag_results$worst_ess_bulk, diag_results$worst_ess_tail),
            divergence_pct = if (total_post_warmup_samples > 0) (sum(num_divergent) / total_post_warmup_samples) * 100 else 0
          )
          
          restart_count <- restart_count + 1
          
          # Calculate new settings for *next* loop
          settings <- calculate_restart_settings(original_warmup, original_adapt_delta, 
                                                 original_max_treedepth, original_n_chains, 
                                                 restart_count, problem_type)
          current_warmup <- settings$warmup
          current_adapt_delta <- settings$adapt_delta
          current_max_treedepth <- settings$max_treedepth
          current_n_chains <- settings$n_chains
          
          cat(sprintf("  Restarting (Attempt %d/%d) with settings: warmup=%d, adapt_delta=%.3f, treedepth=%d, chains=%d\n", 
                      restart_count, MAX_RESTARTS, current_warmup, current_adapt_delta, current_max_treedepth, current_n_chains))
          next  # Go to next restart iteration
          
        } else {
          cat("  MAX_RESTARTS reached. Failing.\n")
          break # Exit restart loop
        }
      }
      
      # 3. If we are here, problem_type == "extension_needed"
      # No fundamental problem - try extending iterations
      while (current_iter < target_total_max && !diag_results$all_pass) {
        
        # --- NEW ESS-based Increment Logic ---
        worst_ess <- min(diag_results$worst_ess_bulk, diag_results$worst_ess_tail, na.rm = TRUE)
        # Use the PASSED-IN target ESS
        target_ess <- diag_thresholds$ess_bulk 
        
        # Default to the small, user-supplied increment
        next_increment <- iter_increment 
        
        if (target_ess >= 800) {
          if (worst_ess < 400) { 
            cat(sprintf("  ESS %d < 400. Applying large increment (%d).\n", floor(worst_ess), initial_post_warmup_chunk))
            next_increment <- initial_post_warmup_chunk 
          } else if (worst_ess < (target_ess / 2)) {
            cat(sprintf("  ESS %d < target/2 (%d). Applying large increment (%d).\n", floor(worst_ess), (target_ess/2), initial_post_warmup_chunk))
            next_increment <- initial_post_warmup_chunk
          }
        }
        # --- End of NEW Logic ---
        
        cat(sprintf("  Diagnostics failed. Adding %d iterations...\n", next_increment))
        target_iter <- min(current_iter + next_increment, target_total_max)
        
        # Checkpoint loop to reach new target
        while (current_iter < target_iter) {
          remaining_iter <- min(checkpoint_interval, target_iter - current_iter)
          
          sampling_result <- continue_sampling(
            stanmodel_arg, data_list, remaining_iter, current_n_chains,
            current_adapt_delta, current_max_treedepth,
            accumulated_samples, step_size, inv_metric
          )
          
          fit <- sampling_result$fit
          new_samples <- sampling_result$new_samples
          current_iter <- current_iter + sampling_result$iterations_added
          
          accumulated_samples <- c(accumulated_samples, list(new_samples))
          new_diagnostics <- fit$sampler_diagnostics()
          accumulated_diagnostics <- c(accumulated_diagnostics, list(new_diagnostics))
          
          save_checkpoint(checkpoint_file, current_iter, accumulated_samples,
                          accumulated_diagnostics, step_size, inv_metric, warmup_done)
          
          cat("Checkpoint saved at iteration", current_iter, "\n")
        }
        
        # Check diagnostics again
        cat(sprintf("Iteration %d complete. Checking diagnostics...\n", current_iter))
        
        all_samples <- do.call(posterior::bind_draws, c(accumulated_samples, along = "iteration"))
        all_diagnostics <- do.call(posterior::bind_draws, c(accumulated_diagnostics, along = "iteration"))
        
        diag_results <- check_diagnostics(all_samples, diag_thresholds)
        
        cat(sprintf("  Worst Rhat: %.3f, Min ESS: %.0f\n",
                    diag_results$worst_rhat,
                    min(c(diag_results$worst_ess_bulk, diag_results$worst_ess_tail), na.rm = TRUE)))
        
        extension_count <- extension_count + 1
      }
      
      # Check final status after extensions
      if (diag_results$all_pass) {
        cat("  All diagnostics passed.\n")
        converged <- TRUE
      } else {
        cat(sprintf("  WARNING: Reached maximum total iterations (%d) without full convergence.\n", target_total_max))
      }
      
      break  # Exit restart loop (either converged or hit max_iter)
    }
    
    # Calculate actual post-warmup iterations
    actual_postwarmup_iter <- current_iter - current_warmup
    
    # Process final results with ACTUAL settings used
    fit_result <- process_sampling_results(
      accumulated_samples, accumulated_diagnostics, current_warmup, actual_postwarmup_iter, current_n_chains,
      current_adapt_delta, current_max_treedepth, model_str, task, n_subs, model_params, output_file
    )
    
    # Add convergence info
    fit_result$convergence_info <- list(
      adaptive_iter_enabled = TRUE,
      converged = converged,
      requested_min_iter = min_iter,
      requested_max_iter = max_iter,
      actual_postwarmup_iter = actual_postwarmup_iter,
      actual_total_iter = current_iter,
      actual_n_chains = current_n_chains,
      extensions_applied = extension_count,
      restart_count = restart_count,
      restart_attempts = if (length(restart_attempts) > 0) restart_attempts else NULL,
      final_worst_rhat = diag_results$worst_rhat,
      final_min_ess_bulk = diag_results$worst_ess_bulk,
      final_min_ess_tail = diag_results$worst_ess_tail
    )
    
  } else {
    # === MCMC Fitting with Checkpointing (Original Behavior) ===
    
    # Initialize or load checkpoint state
    checkpoint_state <- initialize_or_load_checkpoint(checkpoint_file)
    
    # Extract checkpoint state
    current_iter <- checkpoint_state$current_iter
    accumulated_samples <- checkpoint_state$accumulated_samples
    accumulated_diagnostics <- checkpoint_state$accumulated_diagnostics
    step_size <- checkpoint_state$step_size
    inv_metric <- checkpoint_state$inv_metric
    warmup_done <- checkpoint_state$warmup_done
    
    # Total iterations to run
    total_iter <- n_warmup + n_iter
    
    # Run sampling with checkpointing
    while (current_iter < total_iter) {
      # Calculate remaining iterations for this checkpoint interval
      remaining_iter <- min(checkpoint_interval, total_iter - current_iter)
      
      # Initial run or warmup not complete
      if (!warmup_done) {
        sampling_result <- run_initial_sampling(
          stanmodel_arg, data_list, n_warmup, remaining_iter, n_chains, 
          adapt_delta, max_treedepth, init_params
        )
        
        # Extract results
        fit <- sampling_result$fit
        new_samples <- sampling_result$new_samples
        current_iter <- sampling_result$current_iter
        
        # Extract step size and inverse metric after warmup
        step_size <- fit$metadata()$step_size_adaptation
        inv_metric <- fit$inv_metric()
        warmup_done <- TRUE
      } else {
        # Continue sampling post-warmup
        sampling_result <- continue_sampling(
          stanmodel_arg, data_list, remaining_iter, n_chains, adapt_delta, max_treedepth,
          accumulated_samples, step_size, inv_metric
        )
        
        # Extract results
        fit <- sampling_result$fit
        new_samples <- sampling_result$new_samples
        current_iter <- current_iter + sampling_result$iterations_added
      }
      
      # Accumulate new samples and diagnostics
      accumulated_samples <- c(accumulated_samples, list(new_samples))
      new_diagnostics <- fit$sampler_diagnostics()
      accumulated_diagnostics <- c(accumulated_diagnostics, list(new_diagnostics))
      
      # Save checkpoint
      save_checkpoint(
        checkpoint_file, current_iter, accumulated_samples,
        accumulated_diagnostics, step_size, inv_metric, warmup_done
      )
      
      cat("Checkpoint saved at iteration", current_iter, "\n")
    }
    
    # Process final results
    fit_result <- process_sampling_results(
      accumulated_samples, accumulated_diagnostics, n_warmup, n_iter, n_chains,
      adapt_delta, max_treedepth, model_str, task, n_subs, model_params, output_file
    )
    
  }  # End of MCMC branch
  
  # Add additional metadata
  fit_result$cmdstan_version <- cmdstan_version()
  
  # Extract and store subject IDs (REQUIRED)
  if (is.null(data_list$sid)) {
    stop("data_list$sid is required but missing. All models now require subject IDs to be passed.")
  }
  fit_result$subject_list <- data_list$sid
  
  # Store data parameters
  fit_result$data_filt_info = data_filt_list
  
  # Store data quality metrics if available
  if (!is.null(data_list$data_quality)) {
    fit_result$data_quality <- data_list$data_quality
  }
  
  # Save fitted model
  cat("Saving fitted model to:", output_file, "\n")
  saveRDS(fit_result, file = output_file)
  
  # Remove checkpoint file after successful completion
  if (file.exists(checkpoint_file) && file.exists(output_file)) {
    file.remove(checkpoint_file)
  }
  
  return(fit_result)
}

#' Initialize or load checkpoint state
#' @param checkpoint_file Character string with path to checkpoint file
#' @return List with checkpoint state
initialize_or_load_checkpoint <- function(checkpoint_file) {
  if (file.exists(checkpoint_file)) {
    cat("Resuming from checkpoint\n")
    checkpoint <- readRDS(checkpoint_file)
    return(checkpoint)
  } else {
    return(list(
      current_iter = 0,
      accumulated_samples = list(),
      accumulated_diagnostics = list(),
      step_size = NULL,
      inv_metric = NULL,
      warmup_done = FALSE
    ))
  }
}

#' Run initial sampling (including warmup)
#' @param stanmodel_arg Stan model object
#' @param data_list List of data for the model
#' @param n_warmup Integer number of warmup iterations
#' @param remaining_iter Integer number of remaining iterations
#' @param n_chains Integer number of chains
#' @param adapt_delta Numeric adaptation delta parameter
#' @param max_treedepth Integer maximum tree depth
#' @param init_params Optional list of parameter initializations
#' @return List with fit object, new samples, and current iteration
run_initial_sampling <- function(stanmodel_arg, data_list, n_warmup, remaining_iter, n_chains, 
                                 adapt_delta, max_treedepth, init_params) {
  # Ensure iter_to_run is at least n_warmup + 1
  iter_to_run <- max(n_warmup + 1, remaining_iter)
  
  # Adjust init_params if n_chains has changed
  if (!is.null(init_params) && length(init_params) != n_chains) {
    cat("Adjusting init_params list from", length(init_params), "to", n_chains, "chains.\n")
    # Replicate the first init list to match n_chains
    init_params <- replicate(n_chains, init_params[[1]], simplify = FALSE)
  }
  
  # Run sampling
  fit <- stanmodel_arg$sample(
    data = data_list,
    iter_sampling = iter_to_run - n_warmup,
    iter_warmup = n_warmup,
    chains = n_chains,
    parallel_chains = n_chains,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = ceiling(iter_to_run/10),
    init = if (!is.null(init_params)) init_params else NULL
  )
  
  # Return results
  return(list(
    fit = fit,
    new_samples = fit$draws(),
    current_iter = iter_to_run
  ))
}

#' Continue sampling after warmup
#' @param stanmodel_arg Stan model object
#' @param data_list List of data for the model
#' @param remaining_iter Integer number of remaining iterations
#' @param n_chains Integer number of chains
#' @param adapt_delta Numeric adaptation delta parameter
#' @param max_treedepth Integer maximum tree depth
#' @param accumulated_samples List of accumulated samples
#' @param step_size Numeric step size parameter
#' @param inv_metric Inverse metric from warmup
#' @return List with fit object, new samples, and iterations added
continue_sampling <- function(stanmodel_arg, data_list, remaining_iter, n_chains, 
                              adapt_delta, max_treedepth, accumulated_samples, 
                              step_size, inv_metric) {
  # Get the last draws for initialization
  last_draws <- accumulated_samples[[length(accumulated_samples)]]
  last_draws <- last_draws[dim(last_draws)[1],,,drop=FALSE]  # Get the last iteration for all chains
  
  # Check if number of chains in last_draws matches n_chains
  if (dim(last_draws)[2] != n_chains) {
    stop(paste("Mismatch in chains. Accumulated samples have", dim(last_draws)[2],
               "chains, but current run asks for", n_chains))
  }
  
  # Create initialization lists from last draws
  init_list <- lapply(1:n_chains, function(chain_idx) create_init_list(last_draws, chain_idx))
  
  # Run sampling
  fit <- stanmodel_arg$sample(
    data = data_list,
    iter_sampling = remaining_iter,
    iter_warmup = 0,
    chains = n_chains,
    parallel_chains = n_chains,
    adapt_delta = adapt_delta,
    metric = 'dense_e',
    max_treedepth = max_treedepth,
    init = init_list,
    inv_metric = inv_metric,
    adapt_engaged = FALSE,
    step_size = step_size,
    refresh = ceiling(remaining_iter/10)
  )
  
  # Return results
  return(list(
    fit = fit,
    new_samples = fit$draws(),
    iterations_added = dim(fit$draws())[1]
  ))
}

#' Save checkpoint to file
#' @param checkpoint_file Character string with path to checkpoint file
#' @param current_iter Integer current iteration
#' @param accumulated_samples List of accumulated samples
#' @param accumulated_diagnostics List of accumulated diagnostics
#' @param step_size Numeric step size parameter
#' @param inv_metric Inverse metric from warmup
#' @param warmup_done Boolean indicating whether warmup is done
#' @return NULL (saves checkpoint to file)
save_checkpoint <- function(checkpoint_file, current_iter, accumulated_samples,
                            accumulated_diagnostics, step_size, inv_metric, warmup_done) {
  saveRDS(
    list(
      current_iter = current_iter, 
      accumulated_samples = accumulated_samples,
      accumulated_diagnostics = accumulated_diagnostics,
      step_size = step_size, 
      inv_metric = inv_metric,
      warmup_done = warmup_done
    ), 
    file = checkpoint_file
  )
}

#' Process sampling results to create a fit object
#' @param accumulated_samples List of accumulated samples
#' @param accumulated_diagnostics List of accumulated diagnostics
#' @param n_warmup Integer number of warmup iterations
#' @param n_iter Integer number of sampling iterations
#' @param n_chains Integer number of chains
#' @param adapt_delta Numeric adaptation delta parameter
#' @param max_treedepth Integer maximum tree depth
#' @param model_str Character string with model name
#' @param task Character string with task name
#' @param n_subs Integer number of subjects
#' @param model_params Character vector of model parameters
#' @param output_file Character string with output file path
#' @return List with processed sampling results
process_sampling_results <- function(accumulated_samples, accumulated_diagnostics, n_warmup, n_iter, 
                                     n_chains, adapt_delta, max_treedepth, model_str, task,
                                     n_subs, model_params, output_file) {
  # Combine all accumulated samples into a single draws object
  all_samples <- do.call(posterior::bind_draws, c(accumulated_samples, along = "iteration"))
  
  # Combine all accumulated diagnostics
  all_diagnostics <- do.call(posterior::bind_draws, c(accumulated_diagnostics, along = "iteration"))
  
  # Calculate diagnostic summaries
  num_divergent <- colSums(all_diagnostics[, , "divergent__"])
  num_max_treedepth <- colSums(all_diagnostics[, , "treedepth__"] == max_treedepth)
  ebfmi <- sapply(1:n_chains, function(chain) {
    energy <- all_diagnostics[, chain, "energy__"]
    # Handle chains with no variance in energy
    if (var(energy) == 0) return(NA)
    sum(diff(energy)^2) / (length(energy) - 1) / var(energy)
  })
  
  # Calculate summary statistics
  summary_stats <- posterior::summarize_draws(all_samples)
  
  # Calculate parameter diagnostics
  param_names <- dimnames(all_samples)[[3]]
  diagnostics <- calculate_param_diagnostics_batch(all_samples, param_names)
  
  # Print diagnostic warnings
  print_diagnostic_warnings(num_divergent, num_max_treedepth, n_chains, n_iter, n_warmup, max_treedepth)
  
  # Create a fit object with all necessary information
  fit_result = list(
    draws = all_samples,
    sampler_diagnostics = all_diagnostics,
    n_warmup = n_warmup,
    n_iter = n_iter, # This is post-warmup iter
    n_params = dim(all_samples)[3],
    n_chains = n_chains,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    tss = dim(all_samples)[1],
    all_params = param_names,
    list_params = unique(gsub("_subj", "", gsub("\\[.*?\\]", "", param_names))),
    summary_stats = summary_stats,
    diagnostic_summary = list(
      num_divergent = num_divergent,
      num_max_treedepth = num_max_treedepth,
      ebfmi = ebfmi
    ),
    diagnostics = diagnostics,
    model_name = model_str,
    task = task,
    n_runs = length(accumulated_samples),
    filename = output_file
  )
  
  # Extract parameters for diagnostics
  cat("Extracting parameters\n")
  fit_result$params <- extract_params(fit_result$all_params, n_subs = n_subs, main_params_vec = model_params)
  fit_result$params <- unname(fit_result$params)
  fit_result$model_params = model_params
  
  return(fit_result)
}

#' Calculate diagnostics for multiple parameters
#' @param all_samples Array of all samples
#' @param param_names Character vector of parameter names
#' @return Data frame with diagnostics for each parameter
calculate_param_diagnostics_batch <- function(all_samples, param_names) {
  # Apply diagnostic calculation to each parameter
  diagnostics <- lapply(param_names, function(param) {
    param_draws <- all_samples[, , param]
    calculate_param_diagnostics(param_draws)
  })
  
  # Combine the results into a data frame
  diagnostics_df <- do.call(rbind, diagnostics)
  rownames(diagnostics_df) <- param_names
  
  return(diagnostics_df)
}

#' Calculate diagnostics for a single parameter
#' @param param_draws Array of parameter draws
#' @return Vector with diagnostic statistics
calculate_param_diagnostics <- function(param_draws) {
  # Calculate diagnostics
  rhat <- posterior::rhat(param_draws)
  ess_bulk <- posterior::ess_bulk(param_draws)
  ess_tail <- posterior::ess_tail(param_draws)
  
  return(c(rhat = rhat, ess_bulk = ess_bulk, ess_tail = ess_tail))
}

#' Print diagnostic warnings to console
#' @param num_divergent Vector with number of divergent transitions per chain
#' @param num_max_treedepth Vector with number of max treedepth hits per chain
#' @param n_chains Integer number of chains
#' @param n_iter Integer number of sampling iterations
#' @param n_warmup Integer number of warmup iterations
#' @param max_treedepth Integer maximum tree depth
#' @return NULL (prints warnings to console)
print_diagnostic_warnings <- function(num_divergent, num_max_treedepth, 
                                      n_chains, n_iter, n_warmup, max_treedepth) {
  # n_iter is POST-WARMUP iterations
  total_iterations <- n_chains * n_iter 
  total_divergent <- sum(num_divergent)
  total_max_treedepth <- sum(num_max_treedepth)
  
  if (total_iterations == 0) {
    if (total_divergent > 0 || total_max_treedepth > 0) {
      cat("Warning: Diagnostics (divergences/treedepth) occurred, but no post-warmup iterations were saved.\n")
    }
    return()
  }
  
  if (total_divergent > 0) {
    cat(sprintf("\nWarning: %d of %d (%.1f%%) transitions ended with a divergence.\n", 
                total_divergent, total_iterations, 100 * total_divergent / total_iterations))
    cat("See https://mc-stan.org/misc/warnings for details.\n")
  }
  
  if (total_max_treedepth > 0) {
    cat(sprintf("\nWarning: %d of %d (%.1f%%) transitions hit the maximum treedepth limit of %d.\n", 
                total_max_treedepth, total_iterations, 100 * total_max_treedepth / total_iterations, max_treedepth))
    cat("See https://mc-stan.org/misc/warnings for details.\n")
  }
}

#' Validate empirical Bayes fits against hierarchical fits
#' @param hier_fit List with hierarchical fit
#' @param emp_fit List with empirical Bayes fit
#' @param model_params Character vector of model parameters
#' @param subs_df Data frame with subject information
#' @param n_samples Integer number of samples to draw
#' @return List with validation results
validate_empirical_bayes <- function(hier_fit, emp_fit, model_params, subs_df, n_samples = 1000) {
  # Load required helper functions
  source(file.path(here::here(), "scripts", "helper_models.R"))
  
  validation_results <- list()
  
  # Get subjects used in hierarchical fit
  hier_subs <- subs_df$sid[subs_df$set == "hier"]
  
  # Get indices for hierarchical and empirical Bayes fits
  hier_indices <- which(subs_df$set == "hier")
  emp_indices <- which(subs_df$sid %in% hier_subs & subs_df$use == "training")
  
  for (param in model_params) {
    # Extract subject-level parameters for both fits
    hier_params <- hier_fit$draws[,, grep(paste0(param, "_pr\\["), hier_fit$all_params)]
    emp_params <- emp_fit$draws[,, grep(paste0(param, "_pr\\["), emp_fit$all_params)]
    
    # Ensure we're comparing the same subjects
    n_subjects <- length(hier_subs)
    
    # Sample from posterior distributions
    sample_indices <- sample(nrow(hier_params), ceiling(n_samples/2), replace = TRUE)
    
    diffs_pr <- matrix(nrow = ceiling(n_samples/2)*2, ncol = n_subjects)
    diffs <- matrix(nrow = ceiling(n_samples/2)*2, ncol = n_subjects)
    for (i in 1:n_subjects) {
      hier_samples <- hier_params[sample_indices, , paste0(param, "_pr[", hier_indices[i], "]")]
      emp_samples <- emp_params[sample_indices, , paste0(param, "_pr[", emp_indices[i], "]")]
      diffs_pr[, i] <- as.vector(emp_samples - hier_samples)
      
      diffs[, i] <- as.vector(param_xfm(param)(emp_samples) - param_xfm(param)(hier_samples))
    }
    
    # Calculate summary statistics
    validation_results[[paste0(param, "_pr")]] <- list(
      mean_diff = colMeans(diffs_pr),
      median_diff = apply(diffs_pr, 2, median),
      sd_diff = apply(diffs_pr, 2, sd),
      quantiles = t(apply(diffs_pr, 2, quantile, probs = c(0.025, 0.25, 0.75, 0.975)))
    )
    
    # Create a density plot of differences
    plot_data <- data.frame(diff = as.vector(diffs_pr))
    p <- ggplot(plot_data, aes(x = diff)) +
      geom_density(fill = "blue", alpha = 0.5) +
      ggtitle(paste("Difference in", paste0(param, "_pr"), "estimates")) +
      xlab("Empirical Bayes - Hierarchical") +
      theme_minimal()
    
    validation_results[[paste0(param, "_pr")]]$plot <- p
    
    # Add subject IDs to results
    validation_results[[paste0(param, "_pr")]]$subjects <- hier_subs[1:n_subjects]
    
    validation_results[[param]] <- list(
      mean_diff = colMeans(diffs),
      median_diff = apply(diffs, 2, median),
      sd_diff = apply(diffs, 2, sd),
      quantiles = t(apply(diffs, 2, quantile, probs = c(0.025, 0.25, 0.75, 0.975)))
    )
    
    # Create a density plot of differences
    plot_data <- data.frame(diff = as.vector(diffs))
    p <- ggplot(plot_data, aes(x = diff)) +
      geom_density(fill = "blue", alpha = 0.5) +
      ggtitle(paste("Difference in", param, "estimates")) +
      xlab("Empirical Bayes - Hierarchical") +
      theme_minimal()
    
    validation_results[[param]]$plot <- p
    
    # Add subject IDs to results
    validation_results[[param]]$subjects <- hier_subs[1:n_subjects]
  }
  
  return(validation_results)
}



#' Run Pathfinder and process results
#' @param stanmodel_arg Stan model object
#' @param data_list List of data for the model
#' @param num_paths Integer number of pathfinder paths
#' @param draws Integer number of final draws after PSIS
#' @param single_path_draws Integer draws per single path
#' @param psis_resample Boolean whether to use PSIS resampling
#' @param calculate_lp Boolean whether to calculate log probability
#' @param model_str Character string with model name
#' @param task Character string with task name
#' @param n_subs Integer number of subjects
#' @param model_params Character vector of model parameters
#' @param output_file Character string with output file path
#' @param init_params Optional initialization
#' @return List with processed pathfinder results
run_pathfinder_and_process <- function(stanmodel_arg, data_list, num_paths, draws, 
                                       single_path_draws, psis_resample, calculate_lp,
                                       model_str, task, n_subs, model_params, 
                                       output_file, init_params) {
  
  cat("Running Pathfinder with", num_paths, "paths...\n")
  
  # Run Pathfinder
  fit <- stanmodel_arg$pathfinder(
    data = data_list,
    num_paths = num_paths,
    draws = draws,
    single_path_draws = single_path_draws,
    psis_resample = psis_resample,
    calculate_lp = calculate_lp,
    init = if (!is.null(init_params)) init_params else NULL,
    refresh = 100
  )
  
  cat("Pathfinder complete. Processing results...\n")
  
  # Extract draws
  all_samples <- fit$draws()
  
  # Get parameter names
  param_names <- dimnames(all_samples)[[3]]
  
  # Calculate summary statistics
  summary_stats <- posterior::summarize_draws(all_samples)
  
  # Calculate parameter diagnostics
  diagnostics <- calculate_param_diagnostics_batch(all_samples, param_names)
  
  # Extract Pathfinder-specific diagnostics
  pf_diagnostics <- list(
    num_paths = num_paths,
    total_draws = draws,
    method = "pathfinder"
  )
  
  # Extract diagnostic draws (lp__ and lp_approx__)
  diagnostic_vars <- c("lp__", "lp_approx__")
  available_diagnostics <- diagnostic_vars[diagnostic_vars %in% param_names]
  
  if (length(available_diagnostics) > 0) {
    # Create a diagnostics draws object similar to MCMC's sampler_diagnostics
    pf_sampler_diagnostics <- all_samples[, , available_diagnostics, drop = FALSE]
  } else {
    pf_sampler_diagnostics <- NULL
  }
  
  # Extract summary statistics for diagnostics
  if ("lp__" %in% param_names) {
    lp_draws <- all_samples[,,"lp__"]
    pf_diagnostics$lp_mean <- mean(lp_draws)
    pf_diagnostics$lp_sd <- sd(lp_draws)
    pf_diagnostics$lp_min <- min(lp_draws)
    pf_diagnostics$lp_max <- max(lp_draws)
  }
  
  # Extract ELBO statistics (lp_approx__ is the ELBO from Pathfinder)
  if ("lp_approx__" %in% param_names) {
    lp_approx_draws <- all_samples[,,"lp_approx__"]
    pf_diagnostics$elbo_mean <- mean(lp_approx_draws)
    pf_diagnostics$elbo_sd <- sd(lp_approx_draws)
    pf_diagnostics$elbo_min <- min(lp_approx_draws)
    pf_diagnostics$elbo_max <- max(lp_approx_draws)
  }
  
  # Create fit result structure compatible with MCMC output
  fit_result <- list(
    draws = all_samples,
    sampler_diagnostics = pf_sampler_diagnostics,  # Pathfinder diagnostic draws (lp__, lp_approx__)
    n_warmup = 0,  # No warmup in Pathfinder
    n_iter = draws,
    n_params = dim(all_samples)[3],
    n_chains = num_paths,  # Treat paths as "chains" for compatibility
    tss = dim(all_samples)[1],
    all_params = param_names,
    list_params = unique(gsub("_subj", "", gsub("\\[.*?\\]", "", param_names))),
    summary_stats = summary_stats,
    diagnostic_summary = pf_diagnostics,
    diagnostics = diagnostics,
    model_name = model_str,
    task = task,
    fitting_method = "pathfinder",
    filename = output_file
  )
  
  # Extract parameters for diagnostics
  cat("Extracting parameters\n")
  fit_result$params <- extract_params(fit_result$all_params, n_subs = n_subs, main_params_vec = model_params)
  fit_result$params <- unname(fit_result$params)
  fit_result$model_params = model_params
  
  # Print convergence info
  cat("\nPathfinder Diagnostics:\n")
  cat("  Num paths:", num_paths, "\n")
  cat("  Total draws:", draws, "\n")
  if (!is.null(pf_diagnostics$lp_mean)) {
    cat("  Log probability (lp__):\n")
    cat("    Mean:", round(pf_diagnostics$lp_mean, 2), "\n")
    cat("    SD:", round(pf_diagnostics$lp_sd, 2), "\n")
    cat("    Range: [", round(pf_diagnostics$lp_min, 2), ",", round(pf_diagnostics$lp_max, 2), "]\n")
  }
  if (!is.null(pf_diagnostics$elbo_mean)) {
    cat("  ELBO (lp_approx__):\n")
    cat("    Mean:", round(pf_diagnostics$elbo_mean, 2), "\n")
    cat("    SD:", round(pf_diagnostics$elbo_sd, 2), "\n")
    cat("    Range: [", round(pf_diagnostics$elbo_min, 2), ",", round(pf_diagnostics$elbo_max, 2), "]\n")
  }
  
  return(fit_result)
}

create_param_init_list <- function(model_params, no_suffix = NULL, exclude = NULL, default_value = 0) {
  # Remove excluded parameters
  if (!is.null(exclude)) {
    model_params <- model_params[!model_params %in% exclude]
  }
  
  # Create the list
  init_list <- lapply(model_params, function(param) {
    if (param %in% no_suffix) {
      # No suffix for these parameters
      default_value
    } else {
      # Add _pr suffix
      default_value
    }
  })
  
  # Name the list elements
  names(init_list) <- ifelse(model_params %in% no_suffix, 
                             model_params, 
                             paste0(model_params, "_pr"))
  
  return(init_list)
}


create_init_list = function(last_draws, chain_idx) {
  
  # Extract the last draw for this chain
  # last_draws is expected to be [1, n_chains, n_params]
  last_draw_chain <- last_draws[1, chain_idx, ]
  
  # Get all the parameter names
  param_names <- posterior::variables(last_draw_chain)
  if (is.null(param_names)) {
    stop("Last draws are missing parameter names. Cannot create init list.")
  }
  
  # Separate parameters with and without square brackets
  grouped_params <- gsub("\\[.*\\]", "", param_names)  # Remove everything in square brackets
  unique_grouped_params <- unique(grouped_params)      # Get unique base parameter names
  
  # Create the list dynamically
  init_vals <- list()
  
  # Loop through each unique base parameter name
  for (param in unique_grouped_params) {
    # Find all elements corresponding to this parameter
    matching_indices <- which(grouped_params == param)
    
    if (length(matching_indices) > 1) {
      # If there are multiple entries, it's an array/vector: group them
      init_vals[[param]] <- last_draw_chain[matching_indices]
    } else {
      # If there's only one entry, it's a scalar
      init_vals[[param]] <- last_draw_chain[matching_indices]
    }
  }
  
  return(init_vals)
}
