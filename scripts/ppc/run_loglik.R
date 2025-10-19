#!/usr/bin/env Rscript

#' Run Log-likelihood Calculations for Posterior Predictive Checks
#' @description Calculate log-likelihood and information criteria

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(dplyr)
  library(loo)
  library(foreign)  # For read.spss
})

# Define command line options
option_list = list(
  make_option(c("-m", "--model"), type="character", help="Model name (e.g., ev, pvldelta)"),
  make_option(c("-k", "--task"), type="character", help="Task name"),
  make_option(c("-c", "--cohort"), type="character", help="Cohort identifier"),
  make_option(c("-g", "--group"), type="character", default="batch_001", help="Group identifier for file naming"),
  make_option(c("--ses"), type="character", default="00", help="Session identifier"),
  make_option(c("-s", "--sim_file"), type="character", default=NULL, 
              help="Path to simulation file (optional - if not provided, loads data directly)"),
  make_option(c("-f", "--fit_file"), type="character", default=NULL,
              help="Path to fit file (optional - required if sim_file not provided)"),
  make_option(c("-e", "--exclude_file"), type="character", default=NULL, 
              help="File with subject IDs to exclude"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory (optional)"),
  make_option(c("-i", "--ic_method"), type="character", default="loo", 
              help="Information criterion method: loo or waic"),
  make_option(c("-r", "--rt_method"), type="character", default="mark", 
              help="RT handling method: all, remove, force, adaptive, mark"),
  make_option(c("--RTbound_min_ms"), type="numeric", default=50, 
              help="RT lower bound in milliseconds"),
  make_option(c("--RTbound_max_ms"), type="numeric", default=4000, 
              help="RT upper bound in milliseconds"),
  make_option(c("--n_samples"), type="integer", default=100,
              help="Number of posterior samples to use for loglik calculation")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Set default for n_samples if NULL (optparse sometimes doesn't apply defaults)
if (is.null(opt$n_samples)) {
  opt$n_samples <- 100
}

# Check required arguments
if (is.null(opt$model) || is.null(opt$task) || is.null(opt$cohort)) {
  stop("Model name, task name, and cohort are all required.")
}

# Import helper modules - use unified helper entry point
source(file.path(here::here(), "scripts", "ppc", "helpers", "helper_functions_ppc.R"))
source(file.path(here::here(), "scripts", "parameter_recovery", "helper_functions_PR.R"))
source(file.path(here::here(), "scripts", "ppc", "loglik_functions.R"))
source(file.path(here::here(), "scripts", "ppc", "helpers", "helper_ppc_simulation.R"))

# Ensure PPC directories exist
directory_info <- ensure_ppc_dirs(opt$task, opt$cohort, opt$ses)

# Set output directory
output_dir <- if(!is.null(opt$output_dir)) opt$output_dir else get_ppc_loglik_dir(opt$task, opt$cohort, opt$ses)

# Determine data loading mode
use_simulation_file <- !is.null(opt$sim_file)
use_direct_load <- !is.null(opt$fit_file)

if (!use_simulation_file && !use_direct_load) {
  # Try to find simulation file with standard path
  sim_file <- get_ppc_sim_file_path(opt$task, opt$model, opt$group, opt$cohort, opt$ses)
  if (file.exists(sim_file)) {
    message("Found simulation file, using that: ", sim_file)
    use_simulation_file <- TRUE
    opt$sim_file <- sim_file
  } else {
    # Try to find fit file with standard path
    fit_dir <- get_fits_output_dir(opt$task, "fit", opt$cohort, opt$ses)
    group_identifier <- if(opt$group == "hier") "hier" else opt$group
    fit_file <- file.path(fit_dir, 
                         generate_bids_filename(NULL, opt$task, group_identifier, opt$model,
                                               ext = "rds", cohort = opt$cohort, ses = opt$ses,
                                               additional_tags = list("type" = "fit", "desc" = "output")))
    if (file.exists(fit_file)) {
      message("Found fit file, using direct data loading: ", fit_file)
      use_direct_load <- TRUE
      opt$fit_file <- fit_file
    } else {
      stop("Neither simulation file nor fit file provided or found. Specify --sim_file or --fit_file")
    }
  }
}

if (use_simulation_file && use_direct_load) {
  message("Both simulation file and fit file provided. Using simulation file.")
  use_direct_load <- FALSE
}

# Load exclude list if provided
exclude_subjects <- NULL
if (!is.null(opt$exclude_file) && file.exists(opt$exclude_file)) {
  exclude_subjects <- readLines(opt$exclude_file)
  message(paste("Excluding", length(exclude_subjects), "subjects"))
}

# Extract observed data and parameter sets
observed_data_list <- list()
parameter_sets_by_subject <- list()

if (use_simulation_file) {
  # MODE A: Load from simulation file (existing behavior)
  message("Loading observed data and parameters from simulation file...")
  simulation_data <- load_rds_safe(opt$sim_file, "PPC simulation data")
  if (is.null(simulation_data)) {
    stop("Failed to load simulation data from: ", opt$sim_file)
  }
  
  for (subject_id in names(simulation_data)) {
    if (is.null(exclude_subjects) || !subject_id %in% exclude_subjects) {
      subject_sim <- simulation_data[[subject_id]]
      
      # Get observed data
      observed_data_list[[subject_id]] <- subject_sim$observed_data
      
      # Extract parameter sets from simulations
      param_sets <- list()
      n_param_sets <- min(opt$n_samples, length(subject_sim$simulations))
      
      for (i in 1:n_param_sets) {
        param_sets[[i]] <- subject_sim$simulations[[i]]$parameters
      }
      
      # Convert to data frame
      param_df <- do.call(rbind, lapply(param_sets, as.data.frame))
      parameter_sets_by_subject[[subject_id]] <- param_df
    }
  }
  
} else if (use_direct_load) {
  # MODE B: Load observed data and extract parameters from fitted models (new behavior)
  message("Loading observed data directly and extracting parameters from fitted models...")
  
  # Load observed data
  message("Loading observed data...")
  all_data <- load_data(opt$task, opt$cohort, opt$ses)
  
  # Load fitted models
  message("Loading fitted models from: ", opt$fit_file)
  fits <- load_rds_safe(opt$fit_file, "fitted models")
  if (is.null(fits)) {
    stop("Failed to load fitted models from: ", opt$fit_file)
  }
  
  # Get model info
  task_config <- get_task_config(opt$task)
  group_type <- if(opt$group == "hier") "hier" else "sing"
  full_model_name <- paste(opt$task, group_type, opt$model, sep="_")
  model_defaults <- get_model_defaults(opt$task)
  model_params <- model_defaults[[full_model_name]]$params
  
  # Extract parameter sets from posterior
  message("Extracting posterior parameter samples...")
  parameter_sets_by_subject <- extract_posterior_draws(
    fit_file = opt$fit_file,
    model_key = full_model_name,
    model_params = model_params,
    n_samples = opt$n_samples,
    exclude_subjects = exclude_subjects,
    sampling_method = "weighted",  # For loglik, random sampling is fine
    width_control = 0.95
  )
  
  # Extract observed data for each subject
  names(parameter_sets_by_subject) = fits$subject_list
  subject_ids <- names(parameter_sets_by_subject)
  message("Extracting observed data for ", length(subject_ids), " subjects...")
  
  for (subject_id in subject_ids) {
    subject_data <- all_data[all_data$subjID == subject_id, ]
    
    # Convert to list format expected by loglik functions
    if (opt$task == "igt_mod") {
      observed_data_list[[subject_id]] <- list(
        choice = subject_data$choice,
        shown = subject_data$deck,
        deck = subject_data$deck,
        outcome = subject_data$outcome,
        RT = if("RT" %in% names(subject_data)) subject_data$RT else NULL
      )
    } else if (opt$task == "igt") {
      observed_data_list[[subject_id]] <- list(
        choice = subject_data$choice,
        wins = subject_data$wins,
        losses = subject_data$losses,
        RT = if("RT" %in% names(subject_data)) subject_data$RT else NULL
      )
    }
  }
}

# Initialize task and model
SIM_DIR <- file.path(here::here(), "scripts", "simulation")

# Source required base files first
source_required_files(SIM_DIR)

# Initialize task and model
task <- initialize_task(opt$task, SIM_DIR)
model <- initialize_model(opt$model, opt$task, task, SIM_DIR)

# Build task_params from command line RT bounds
task_params <- list(
  RTbound_min = opt$RTbound_min_ms / 1000,  # Convert ms to seconds
  RTbound_max = opt$RTbound_max_ms / 1000
)

message("Using RT bounds: [", task_params$RTbound_min, ", ", 
        task_params$RTbound_max, "] seconds")
message("Using RT method: ", opt$rt_method)

# ONLY preprocess RTs if model has SSM component
has_ssm <- model$model_type %in% c("SSM", "RL_SSM")

if (has_ssm) {
  # CRITICAL: Preprocess observed data with rt_method to match how model was fitted
  message("Model type: ", model$model_type, " - preprocessing observed data with rt_method=", opt$rt_method, "...")
  for (subject_id in names(observed_data_list)) {
    obs_data <- observed_data_list[[subject_id]]
    
    # Check if RT data exists
    if ("RT" %in% names(obs_data)) {
      # Convert to data frame format expected by preprocess_data
      if (opt$task == "igt"){
        obs_df <- data.frame(
          subjID = subject_id,
          RT = obs_data$RT,
          choice = obs_data$choice,
          wins = if ("wins" %in% names(obs_data)) obs_data$wins else NULL,
          losses = if ("losses" %in% names(obs_data)) obs_data$losses else NULL
        )
      } else if (opt$task == "igt_mod"){
        obs_df <- data.frame(
          subjID = subject_id,
          RT = obs_data$RT,
          choice = obs_data$choice,
          shown = if ("shown" %in% names(obs_data)) obs_data$shown else NULL,
          outcome = if ("outcome" %in% names(obs_data)) obs_data$outcome else NULL,
        )
      }
      
      
      # Apply RT preprocessing
      processed_df <- preprocess_data(
        data = obs_df,
        task = opt$task,
        RTbound_min_ms = opt$RTbound_min_ms,
        RTbound_max_ms = opt$RTbound_max_ms,
        rt_method = opt$rt_method,
        return_dropped_indices = FALSE
      )
      
      # Update the observed data with preprocessed RTs
      observed_data_list[[subject_id]]$RT <- processed_df$RT
      
      # If rt_method="remove", we also need to filter other columns
      if (opt$rt_method == "remove" && nrow(processed_df) < nrow(obs_df)) {
        message("  Subject ", subject_id, ": removed ", nrow(obs_df) - nrow(processed_df), 
                " trials with invalid RTs")
        # Filter all columns to match the kept trials
        for (col in names(observed_data_list[[subject_id]])) {
          if (length(observed_data_list[[subject_id]][[col]]) == nrow(obs_df)) {
            observed_data_list[[subject_id]][[col]] <- observed_data_list[[subject_id]][[col]][1:nrow(processed_df)]
          }
        }
      }
    }
  }
} else {
  message("Model type: ", model$model_type, " - skipping RT preprocessing (pure RL model)")
}

# Calculate log-likelihood
message("Calculating log-likelihood...")
loglik_results <- calculate_loglik_multiple_subjects(
  task_name = opt$task,
  model = model,
  subject_data_list = observed_data_list,
  parameter_sets_by_subject = parameter_sets_by_subject,
  task_params = task_params
)

# Calculate information criteria
message("Calculating information criteria...")
ic_results <- compute_information_criteria(loglik_results, opt$ic_method)
combined_ic <- combine_information_criteria(ic_results)

# Add to results
loglik_results$ic_results <- ic_results
loglik_results$combined_ic <- combined_ic
loglik_results$model <- opt$model
loglik_results$group <- opt$group
loglik_results$task <- opt$task
loglik_results$cohort <- opt$cohort
loglik_results$session <- opt$ses

# Save log-likelihood results using helper function
loglik_file <- save_loglik_results(loglik_results, opt$task, opt$group, opt$model, opt$cohort, output_dir, opt$ses)

message("Log-likelihood calculation complete.")
message("Results saved to: ", loglik_file)

# Print basic summary
message("\nLog-likelihood Summary:")
total_loglik <- sum(sapply(loglik_results[names(parameter_sets_by_subject)], function(x) sum(x$total_loglik, na.rm = TRUE)), na.rm = TRUE)
total_obs <- sum(sapply(loglik_results[names(parameter_sets_by_subject)], function(x) x$n_trials * x$n_param_sets), na.rm = TRUE)
mean_loglik <- total_loglik / total_obs

message("Total log-likelihood: ", round(total_loglik, 2))
message("Total observations: ", total_obs)
message("Mean log-likelihood: ", round(mean_loglik, 4))

if (opt$ic_method == "loo") {
  if (!is.null(combined_ic)) {
    loo_estimate <- combined_ic$estimates["looic", "Estimate"]
    loo_se <- combined_ic$estimates["looic", "SE"]
    message("LOOIC: ", round(loo_estimate, 2), " (SE ", round(loo_se, 2), ")")
  } else {
    message("LOOIC: Unable to calculate")
  }
} else if (opt$ic_method == "waic") {
  if (!is.null(combined_ic)) {
    waic_estimate <- combined_ic$estimates["waic", "Estimate"]
    waic_se <- combined_ic$estimates["waic", "SE"]
    message("WAIC: ", round(waic_estimate, 2), " (SE ", round(waic_se, 2), ")")
  } else {
    message("WAIC: Unable to calculate")
  }
}
