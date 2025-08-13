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
#' @param n_warmup Integer number of warmup iterations
#' @param n_iter Integer number of sampling iterations
#' @param n_chains Integer number of chains
#' @param adapt_delta Numeric adaptation delta parameter
#' @param max_treedepth Integer maximum tree depth
#' @param model_params Character vector of model parameters
#' @param dry_run Boolean to perform a dry run without actual sampling
#' @param checkpoint_interval Integer number of iterations between checkpoints
#' @param output_dir Optional character string with output directory path
#' @param emp_bayes Boolean indicating whether this is an empirical Bayes model
#' @param informative_priors Optional list of informative priors
#' @param subid Optional character string with subject ID
#' @param index Optional integer subject index
#' @param init_params Optional list of parameter initializations
#' @return List with fitted model results or NULL for dry run
fit_and_save_model <- function(task, cohort, ses, group_type, model_name, model_type, data_list, 
                              n_subs, n_trials, n_warmup, n_iter, n_chains, adapt_delta, max_treedepth, 
                              model_params, dry_run = FALSE, checkpoint_interval = 1000, 
                              output_dir = NULL, emp_bayes = FALSE, informative_priors = NULL,
                              subid = NULL, subject_list = NULL, index = NULL, init_params = NULL, cohort_sub_dir = TRUE) {
  
  # Create the model string and get the model path
  model_str <- paste(task, group_type, model_name, sep="_")
  model_path <- get_model_file_path(task, group_type, model_name, model_type)
  
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
        stop("Number of initialization lists (", length(init_params), 
             ") must match number of chains (", n_chains, ")")
      }
    } else {
      # If single list, replicate it for each chain
      init_params <- replicate(n_chains, init_params, simplify = FALSE)
    }
  }
  
  # Determine output directory and file paths
  if (is.null(output_dir)) {
    if (emp_bayes) {
      output_dir <- get_empbayes_dir(task)
    } else {
      output_dir <- get_rds_dir(task, model_type)
    }
  }
  
  # Ensure output directory exists
  ensure_dir_exists(output_dir)
  
  # Generate output and checkpoint file paths
  output_file <- get_output_file_path(task, cohort, group_type, model_name, model_type, 
                                     emp_bayes, subid, index, ses = ses, 
                                     output_dir = output_dir, cohort_sub_dir)
  output_file <- fix_duplicate_cohort_path(output_file, cohort)
  checkpoint_file <- paste0(tools::file_path_sans_ext(output_file), "_checkpoint.rds")
  
  # Dry run only prints information without fitting
  if (dry_run) {
    cat("Dry run for model:", model_str, "\n")
    cat("Data list:\n")
    print(str(data_list))
    cat("Parameters:\n")
    cat("n_subs =", n_subs, "\n")
    cat("n_warmup =", n_warmup, "\n")
    cat("n_iter =", n_iter, "\n")
    cat("n_chains =", n_chains, "\n")
    cat("adapt_delta =", adapt_delta, "\n")
    cat("max_treedepth =", max_treedepth, "\n")
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
  
  # Add informative priors to data_list if provided
  if (!is.null(informative_priors)) {
    data_list$pr_mu <- unname(sapply(model_params, function(param) informative_priors[[param]][1]))
    data_list$pr_sigma <- unname(sapply(model_params, function(param) informative_priors[[param]][2]))
  }
  
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
  
  # Add additional metadata
  fit_result$cmdstan_version <- cmdstan_version()
  
  # Add subid and index to the fit object if provided
  if (!is.null(subid)) fit_result$subid <- subid
  if (!is.null(index)) fit_result$index <- index
  if (!is.null(subject_list)) fit_result$subject_list <- subject_list
  
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
  last_draws <- last_draws[dim(last_draws)[1],,]  # Get the last iteration for all chains
  
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
    n_iter = n_iter,
    n_params = dim(all_samples)[3],
    n_chains = n_chains,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    tss = dim(all_samples)[1],
    all_params = param_names,
    list_params = unique(gsub("\\[.*?\\]", "", param_names)),
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
  total_iterations <- n_chains * (n_iter - n_warmup)
  total_divergent <- sum(num_divergent)
  total_max_treedepth <- sum(num_max_treedepth)
  
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
  last_draw_chain <- last_draws[1, chain_idx, ]
  
  # Get all the parameter names
  param_names <- dimnames(last_draw_chain)[[3]]
  
  # Separate parameters with and without square brackets
  grouped_params <- gsub("\\[.*\\]", "", param_names)  # Remove everything in square brackets
  unique_grouped_params <- unique(grouped_params)      # Get unique base parameter names
  
  # Create the list dynamically
  init_vals <- list()
  
  # Loop through each unique base parameter name
  for (param in unique_grouped_params) {
    # Find all elements corresponding to this parameter
    matching_indices <- grep(paste0("^", param), param_names)
    
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

