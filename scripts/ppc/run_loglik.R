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
  make_option(c("--ses"), type="character", default=NULL, help="Session identifier"),
  make_option(c("-s", "--sim_file"), type="character", default=NULL, 
              help="Path to simulation file (optional, built from task, group, model if not provided)"),
  make_option(c("-e", "--exclude_file"), type="character", default=NULL, 
              help="File with subject IDs to exclude"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory (optional)"),
  make_option(c("-i", "--ic_method"), type="character", default="loo", 
              help="Information criterion method: loo or waic"),
  make_option(c("-r", "--rt_method"), type="character", default="remove", 
              help="RT handling method: all, remove, force, adaptive"),
  make_option(c("--RTbound_min_ms"), type="numeric", default=50, 
              help="RT lower bound in milliseconds"),
  make_option(c("--RTbound_max_ms"), type="numeric", default=4000, 
              help="RT upper bound in milliseconds")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$model) || is.null(opt$task) || is.null(opt$cohort)) {
  stop("Model name, task name, and cohort are all required.")
}

# Import helper modules
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "parameter_recovery", "helper_functions_PR.R"))
source(file.path(here::here(), "scripts", "ppc", "helpers", "helper_ppc_dirs.R"))
source(file.path(here::here(), "scripts", "ppc", "helpers", "task_config.R"))
source(file.path(here::here(), "scripts", "ppc", "loglik_functions.R"))
source(file.path(here::here(), "scripts", "ppc", "helpers", "helper_ppc_simulation.R"))

# Ensure PPC directories exist
directory_info <- ensure_ppc_dirs(opt$task, opt$cohort, opt$ses)

# Set output directory
output_dir <- if(!is.null(opt$output_dir)) opt$output_dir else get_ppc_loglik_dir(opt$task, opt$cohort, opt$ses)

# Build simulation file path if not provided
if (is.null(opt$sim_file)) {
  # Use helper function to get the PPC simulation file path
  sim_file <- get_ppc_sim_file_path(opt$task, opt$model, opt$group, opt$cohort, opt$ses)
  message("Using simulation file: ", sim_file)
} else {
  sim_file <- opt$sim_file
}

# Check if simulation file exists
if (!file.exists(sim_file)) {
  stop("Simulation file not found: ", sim_file)
}

# Load simulation data
message("Loading simulation data...")
simulation_data <- readRDS(sim_file)

# Load exclude list if provided
exclude_subjects <- NULL
if (!is.null(opt$exclude_file) && file.exists(opt$exclude_file)) {
  exclude_subjects <- readLines(opt$exclude_file)
  message(paste("Excluding", length(exclude_subjects), "subjects"))
}

# Extract observed data and parameter sets from simulation results
message("Extracting data from simulation results...")
observed_data_list <- list()
parameter_sets_by_subject <- list()

for (subject_id in names(simulation_data)) {
  if (is.null(exclude_subjects) || !subject_id %in% exclude_subjects) {
    subject_sim <- simulation_data[[subject_id]]
    
    # Get observed data
    observed_data_list[[subject_id]] <- subject_sim$observed_data
    
    # Extract parameter sets from first few simulations
    param_sets <- list()
    n_param_sets <- min(100, length(subject_sim$simulations))  # Use up to 100 parameter sets
    
    for (i in 1:n_param_sets) {
      param_sets[[i]] <- subject_sim$simulations[[i]]$parameters
    }
    
    # Convert to data frame
    param_df <- do.call(rbind, lapply(param_sets, as.data.frame))
    parameter_sets_by_subject[[subject_id]] <- param_df
  }
}

# Initialize task and model
PR_DIR <- file.path(here::here(), "scripts", "parameter_recovery")

# Source required base files first
source_required_files(PR_DIR)

# Initialize task and model
task <- initialize_task(opt$task, PR_DIR)
model <- initialize_model(toupper(opt$model), tolower(opt$task), task, PR_DIR)

# Build task_params from command line RT bounds
task_params <- list(
  RTbound_min = opt$RTbound_min_ms / 1000,  # Convert ms to seconds
  RTbound_max = opt$RTbound_max_ms / 1000
)

message("Using RT bounds: [", task_params$RTbound_min, ", ", 
        task_params$RTbound_max, "] seconds")

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
