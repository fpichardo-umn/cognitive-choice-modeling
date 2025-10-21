#!/usr/bin/env Rscript

#' Run Simulation Stage of PPC Pipeline
#' @description Generate simulated data from model posterior distributions

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(dplyr)
  library(parallel)
  library(cmdstanr)
  library(posterior)
})

# Define command line options
option_list = list(
  make_option(c("-m", "--model"), type="character", default=NULL, help="Model name"),
  make_option(c("-k", "--task"), type="character", default=NULL, help="Task name"),
  make_option(c("-c", "--cohort"), type="character", default=NULL, help="Cohort identifier"),
  make_option(c("-g", "--group"), type="character", default="sing", help="Fit type: sing (individual) or hier (hierarchical)"),
  make_option(c("--group_name"), type="character", default="batch_001", help="Batch identifier (only used for individual fits)"),
  make_option(c("--ses"), type="character", default=NULL, help="Session identifier"),
  make_option(c("-f", "--fit_file"), type="character", default=NULL, 
              help="Path to fit file (optional, built from task, group, model if not provided)"),
  make_option(c("-n", "--n_sims"), type="integer", default=100, 
              help="Number of simulations per subject"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory (optional)"),
  make_option(c("-p", "--parallel"), action="store_true", default=FALSE, 
              help="Use parallel processing"),
  make_option(c("--n_cores"), type="integer", default=2, 
              help="Number of cores for parallel processing"),
  make_option(c("--sampling"), type="character", default="weighted", 
              help="Sampling method for posterior draws: random, width, or weighted"),
  make_option(c("--width_control"), type="numeric", default=0.95, 
              help="Width control parameter for width sampling method (0-1)"),
  make_option(c("--rt_method"), type="character", default="remove", 
              help="RT handling method: all, remove, force, or adaptive"),
  make_option(c("--RTbound_min_ms"), type="numeric", default=50, 
              help="RT lower bound in milliseconds"),
  make_option(c("--RTbound_max_ms"), type="numeric", default=4000, 
              help="RT upper bound in milliseconds"),
  make_option(c("--exclude_file"), type="character", default=NULL, 
              help="Path to file with subject IDs to exclude")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$model) || is.null(opt$task) || is.null(opt$cohort)) {
  stop("Model name, task name, and cohort are all required.")
}

# Validate group parameter
if (!opt$group %in% c("sing", "hier")) {
  stop("Invalid --group parameter. Must be 'sing' (individual) or 'hier' (hierarchical).")
}

# Set up environment - use unified helper entry point
script_dir <- file.path(here::here(), "scripts")
source(file.path(script_dir, "ppc", "helpers", "helper_functions_ppc.R"))
source(file.path(script_dir, "ppc", "simulation_functions.R"))

# Set output directory
if (is.null(opt$output_dir)) {
  output_dir <- get_ppc_sim_dir(opt$task, opt$cohort, opt$ses)
} else {
  output_dir <- opt$output_dir
}

# Load subject exclusion list if provided
exclude_subjects <- NULL
if (!is.null(opt$exclude_file)) {
  if (file.exists(opt$exclude_file)) {
    exclude_subjects <- readLines(opt$exclude_file)
    message(paste("Excluding", length(exclude_subjects), "subjects from", opt$exclude_file))
  } else {
    warning(paste("Exclusion file not found:", opt$exclude_file))
  }
}

# Construct fit file path if not provided
if (is.null(opt$fit_file)) {
  fit_file = file.path(here::here(), "Outputs", opt$task,
            "fits", "fit", opt$cohort, paste0("ses-", opt$ses), 
            generate_bids_filename(
              prefix = NULL,
              task = opt$task,
              group = if(opt$group == "sing") opt$group_name else opt$group,
              model = opt$model,
              cohort = opt$cohort,
              ses = opt$ses,
              additional_tags = list(
                "type" = "fit",
                "desc" = "output"
              ),
              ext = "rds"
            ))
  
  if (!file.exists(fit_file)) {
    stop("Fit file not found: ", fit_file, ". Use --fit_file to specify path.")
  }
} else {
  fit_file <- opt$fit_file
  if (!file.exists(fit_file)) {
    stop("Specified fit file not found: ", fit_file)
  }
}

message("Using fit file: ", fit_file)

# Load fits with error handling
fits <- load_rds_safe(fit_file, "fitted model")
if (is.null(fits)) {
  stop("Failed to load fitted model from: ", fit_file)
}

# Get task configuration
task_config <- get_task_config(opt$task)
model_type <- task_config$type

# Extract model parameters based on task type
model_name <- opt$model
task <- opt$task
group_type <- opt$group  # Use the group parameter (sing or hier)
full_model_name = paste(task, group_type, model_name, sep="_")

model_defaults <- get_model_defaults(opt$task)
if (!full_model_name %in% names(model_defaults)) {
  cat(full_model_name,"\n")
  cat(names(model_defaults))
  stop("Unrecognized model. Please check the model name.")
}

data_to_extract <- if (!is.null(opt$data)) strsplit(opt$data, ",")[[1]] else model_defaults[[full_model_name]]$data
model_params <- if (!is.null(opt$params)) strsplit(opt$params, ",")[[1]] else model_defaults[[full_model_name]]$params

cat("Preparing data for", full_model_name, "\n")


# Load data using the same method as fit scripts
all_data <- load_data(opt$task, opt$cohort, opt$ses)
message("Loaded data with ", nrow(all_data), " trials for ", length(unique(all_data$subjID)), " subjects")

# Extract posterior draws with density-aware sampling
parameter_sets_by_subject <- extract_posterior_draws(
  fit_file = fit_file,
  model_key = full_model_name,  # Pass full model name so hierarchical detection works
  model_params = model_params,
  n_samples = opt$n_sims,
  exclude_subjects = exclude_subjects,
  sampling_method = opt$sampling,
  width_control = opt$width_control
)

message("Loaded fitted model with ", length(parameter_sets_by_subject), " subjects")

if (opt$group == "hier"){
  subject_ids = unique(all_data$subjID)[1:length(names(parameter_sets_by_subject))]
  names(parameter_sets_by_subject) = subject_ids
  
  # Filter data for subjects that have fits
  subject_data <- all_data[all_data$subjID %in% subject_ids, ]
} else {
  # Get subject IDs from the extraction results
  subject_ids <- names(parameter_sets_by_subject)
  
  # Filter data for subjects that have fits
  subject_data <- all_data[all_data$subjID %in% subject_ids, ]
}

message("Using data for ", length(unique(subject_data$subjID)), " subjects with fitted models")

message("Running simulations for ", length(parameter_sets_by_subject), " subjects")

# Build task_params from command line RT bounds
task_params <- list(
  RTbound_min = opt$RTbound_min_ms / 1000,  # Convert ms to seconds
  RTbound_max = opt$RTbound_max_ms / 1000
)

message("Using RT bounds: [", task_params$RTbound_min, ", ", 
        task_params$RTbound_max, "] seconds")

# Generate simulation data
simulation_results <- generate_simulation_data(
  task_name = opt$task,
  model_name = opt$model,
  subject_data = subject_data,
  parameter_sets_by_subject = parameter_sets_by_subject,
  task_params = task_params
)

# Save results with error handling and verification
# Output file naming should reflect the source of parameters:
# - For hierarchical fits: use "hier" (even though we simulate individuals)
# - For individual fits: use group_name (batch identifier)
group_for_output_files <- if(opt$group == "hier") "hier" else opt$group_name
output_file <- get_ppc_sim_file_path(opt$task, opt$model, group_for_output_files, opt$cohort, opt$ses)

success <- save_rds_safe(simulation_results, output_file, "PPC simulation results")
if (!success) {
  stop("Failed to save simulation results to: ", output_file)
}

message("Simulation complete.")
