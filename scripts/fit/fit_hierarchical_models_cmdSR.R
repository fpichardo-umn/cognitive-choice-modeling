#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(cmdstanr)
  library(posterior)
  library(foreign)
  library(dplyr)
  library(tidyr)
})

# Parse command line arguments
option_list = list(
  make_option(c("-m", "--model"), type="character", default=NULL, help="Model name"),
  make_option(c("-t", "--type"), type="character", default="fit", help="Model type (fit, postpc, prepc)"),
  make_option(c("-k", "--task"), type="character", default=NULL, help="Task name"),
  make_option(c("-g", "--group"), type="character", default=NULL, help="Group type (sing, group, group_hier)"),
  make_option(c("-s", "--source"), type="character", default=NULL, help="Data source"),
  make_option(c("--ses"), type="character", default=NULL, help="Session identifier (optional)"),
  make_option(c("--model_status"), type="character", default=NULL,
              help="Model status (canonical/experimental/working). Default: auto-detect"),
  make_option(c("-d", "--data"), type="character", default=NULL, 
              help="Comma-separated list of data to extract"),
  make_option(c("-p", "--params"), type="character", default=NULL, 
              help="Comma-separated list of model parameters"),
  make_option(c("--n_subs"), type="integer", default=1000, help="Number of subjects for hierarchical model"),
  make_option(c("-l", "--subs_file"), type="character", default="subject_ids_complete_valid.txt", 
              help="Subs list file [Data/raw/COHORT/ses-SES/] (default: subject_ids_complete_valid.txt)"),
  make_option(c("--n_trials"), type="integer", default=100, help="Number of trials"),
  make_option(c("--RTbound_min_ms"), type="integer", default=50, help="RT min bound in milliseconds"),
  make_option(c("--RTbound_max_ms"), type="integer", default=4000, help="RT max bound in milliseconds"),
  make_option(c("--rt_method"), type="character", default="mark", help="RT method"),
  make_option(c("--n_warmup"), type="integer", default=3000, help="Number of warmup iterations"),
  make_option(c("--n_iter"), type="integer", default=15000, help="Number of iterations"),
  make_option(c("--n_chains"), type="integer", default=4, help="Number of chains"),
  make_option(c("--adapt_delta"), type="double", default=0.95, help="Adapt delta"),
  make_option(c("--max_treedepth"), type="integer", default=12, help="Max tree depth"),
  make_option(c("--seed"), type="integer", default=NULL, help="Set seed. Default: random"),
  make_option(c("--dry_run"), action="store_true", default=FALSE, help="Perform a dry run"),
  make_option(c("--check_iter"), type="integer", default=1000, help="Iteration interval for checkpoint runs. Default: 1000"),
  make_option(c("--init"), action="store_true", default=FALSE, help="Initialize values to 0"),
  make_option(c("--min_valid_rt_pct"), type="double", default=0.7, help="Minimum percent valid RT"),
  make_option(c("--min_iter"), type="integer", default=NULL, 
              help="Minimum post-warmup iterations for adaptive fitting (enables adaptive mode)"),
  make_option(c("--max_iter"), type="integer", default=NULL,
              help="Maximum post-warmup iterations for adaptive fitting"),
  make_option(c("--iter_increment"), type="integer", default=1000,
              help="Iterations to add at each diagnostic check (default: 1000)"),
  make_option(c("--target_rhat"), type="double", default=1.01,
              help="Target Rhat threshold (default: 1.01)"),
  make_option(c("--target_ess_bulk"), type="integer", default=400,
              help="Target ESS bulk threshold (default: 400)"),
  make_option(c("--target_ess_tail"), type="integer", default=400,
              help="Target ESS tail threshold (default: 400)"),
  make_option(c("--disable_adaptive_iter"), action="store_true", default=FALSE,
              help="Disable adaptive iteration feature (use fixed n_iter)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Options used:\n")
dput(opt)

#set_cmdstan_path("~/stan/cmdstan-2.36.0/")

# Load helper functions for directory structure
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))

# Set up directories using helper functions
PROJ_DIR <- get_proj_dir()
SCRIPT_DIR <- file.path(PROJ_DIR, "scripts")
SAFE_DATA_DIR <- get_safe_data_dir()

# Get task-specific directories
if (is.null(opt$task)) {
  stop("Task name is required using the -k option.")
}

# Set random seed for reproducibility
# If seed not provided, generate one
if (is.null(opt$seed)) {
  opt$seed <- sample.int(.Machine$integer.max, 1)
  message(sprintf("No seed provided. Using random seed: %d", opt$seed))
}
set.seed(opt$seed)

# Get model defaults using helper function
model_defaults <- get_model_defaults(opt$task)

# Main execution - validate hierarchical-specific arguments
if (is.null(opt$model) || is.null(opt$task) || is.null(opt$group) || is.null(opt$source)) {
  stop("Please specify a model, task, group type, and source using the -m, -k, -g, and -s options.")
}

# Validate that this is actually a hierarchical model request
if (!opt$group %in% c("group", "hier")) {
  stop("This script is for hierarchical models. Group type must be 'group' or 'group_hier'. For single subjects, use fit_single_model.R")
}

model_name <- opt$model
task <- opt$task
group_type <- opt$group
full_model_name = paste(task, group_type, model_name, sep="_")

if (!full_model_name %in% names(model_defaults)) {
  cat("Requested model:", full_model_name,"\n")
  cat("Available models:", names(model_defaults), "\n")
  stop("Unrecognized hierarchical model. Please check the model name and group type.")
}

# Extract data and parameter specifications
data_to_extract <- if (!is.null(opt$data)) strsplit(opt$data, ",")[[1]] else model_defaults[[full_model_name]]$data
model_params <- if (!is.null(opt$params)) strsplit(opt$params, ",")[[1]] else model_defaults[[full_model_name]]$params
non_pr_params <- if (opt$init) model_defaults[[full_model_name]]$non_pr_params else NULL
exclude_params <- if (opt$init) model_defaults[[full_model_name]]$exclude_params else NULL

cat("Preparing hierarchical data for", full_model_name, "\n")
cat("Target sample size:", opt$n_subs, "subjects\n")

# Get task-specific directories
if (!is.null(opt$task)) {
  # Get subjects directory with source
  SUBS_DIR <- file.path(SAFE_DATA_DIR, opt$source)
  if (!is.null(opt$ses)) {
    SUBS_DIR <- file.path(SUBS_DIR, paste0("ses-", opt$ses))
  }
  SUBS_LIST_FILE <- file.path(SUBS_DIR, opt$subs_file)
} else {
  stop("Task name is required using the -k option.")
}

# Check subjects list file
if (!file.exists(SUBS_LIST_FILE)) {
  stop(sprintf("Subjects list file not found: %s", SUBS_LIST_FILE))
}
subject_ids <- readLines(SUBS_LIST_FILE)

# Load and prepare hierarchical data
if (!opt$dry_run) {
  cat("Loading data for hierarchical model fitting...\n")
  
  # Load full dataset using helper function
  all_data <- load_data(opt$task, opt$source, opt$ses)
  
  # Filter based on list
  all_data = all_data[all_data$subjID %in% subject_ids,]
  
  # For hierarchical models, we work with multiple subjects
  # Check if we have enough subjects
  unique_subjects <- unique(all_data$subjID)
  n_available_subs <- length(unique_subjects)
  
  cat("Available subjects in dataset:", n_available_subs, "\n")
  
  if (n_available_subs < opt$n_subs) {
    cat("Warning: Requested", opt$n_subs, "subjects but only", n_available_subs, "available. Using all available subjects.\n")
    selected_n_subs <- n_available_subs
  } else {
    selected_n_subs <- opt$n_subs
  }
  
  # For hierarchical models, we typically sample subjects or use all available
  # Here we'll use the first n_subs subjects (could be randomized if needed)
  selected_subjects <- as.character(unique_subjects[1:selected_n_subs])
  hierarchical_data <- all_data[all_data$subjID %in% selected_subjects, ]
  
  cat("Selected", selected_n_subs, "subjects for hierarchical modeling\n")
  
  # Extract data for hierarchical model structure
  data_list <- extract_sample_data(hierarchical_data, data_to_extract, 
                                   task = opt$task,
                                   n_subs = selected_n_subs, 
                                   n_trials = opt$n_trials, 
                                   RTbound_min_ms = opt$RTbound_min_ms, 
                                   RTbound_max_ms = opt$RTbound_max_ms,
                                   RTbound_reject_min_ms = opt$RTbound_min_ms + 20, 
                                   RTbound_reject_max_ms = opt$RTbound_max_ms, 
                                   rt_method = opt$rt_method, 
                                   minrt_ep_ms = 0, min_valid_rt_pct = opt$min_valid_rt_pct)
  
  # Collect data filtering info
  data_filt = c(
    n_trials = opt$n_trials, 
    RTbound_min_ms = opt$RTbound_min_ms, 
    RTbound_max_ms = opt$RTbound_max_ms,
    RTbound_reject_min_ms = opt$RTbound_min_ms + 20, 
    RTbound_reject_max_ms = opt$RTbound_max_ms, 
    rt_method = opt$rt_method, 
    minrt_ep_ms = 0,
    min_valid_rt_pct = opt$min_valid_rt_pct
  )
  
  # Handle entropy data if present (convert matrices to vectors)
  if (sum(grepl("entropy", names(data_list))) > 0){
    cat("Processing entropy measures for hierarchical model...\n")
    if ("shannon_entropy" %in% names(data_list)) {
      data_list$shannon_entropy <- as.vector(data_list$shannon_entropy)
    }
    if ("transition_entropy" %in% names(data_list)) {
      data_list$transition_entropy <- as.vector(data_list$transition_entropy)
    }
    if ("ngram_entropy" %in% names(data_list)) {
      data_list$ngram_entropy <- as.vector(data_list$ngram_entropy)
    }
    if ("conditional_entropy" %in% names(data_list)) {
      data_list$conditional_entropy <- as.vector(data_list$conditional_entropy)
    }
  }
  
} else {
  cat("Dry run: Data would be loaded using load_data(", opt$task, ", ", opt$source, ", ", opt$ses, ")\n")
  cat("Dry run: Would select", opt$n_subs, "subjects for hierarchical modeling\n")
  cat("Dry run: Data to extract:", paste(data_to_extract, collapse=", "), "\n")
  data_list <- NULL
  selected_n_subs <- opt$n_subs
}

# Generate initialization values for hierarchical model
if (opt$init) {
  cat("Creating parameter initialization for hierarchical model...\n")
  model_init_vals = create_param_init_list(model_params, 
                                           no_suffix = non_pr_params, 
                                           exclude = exclude_params)
} else {
  model_init_vals = NULL
}

# Set up output directory for hierarchical models
output_dir <- get_fits_output_dir(opt$task, opt$type, opt$source, opt$ses)

# Construct diagnostic thresholds
diag_thresholds <- list(
  rhat = opt$target_rhat,
  ess_bulk = opt$target_ess_bulk,
  ess_tail = opt$target_ess_tail
)

# Fit hierarchical model
cat("Fitting hierarchical model:", full_model_name, "\n")
fit <- fit_and_save_model(task, opt$source, opt$ses, group_type, model_name,
                          opt$type, data_list, 
                          n_subs = selected_n_subs, 
                          n_trials = opt$n_trials, 
                          n_warmup = opt$n_warmup, 
                          n_iter = opt$n_iter, 
                          n_chains = opt$n_chains,
                          adapt_delta = opt$adapt_delta, 
                          max_treedepth = opt$max_treedepth, 
                          model_params = model_params, 
                          dry_run = opt$dry_run, 
                          checkpoint_interval = opt$check_iter, 
                          output_dir = output_dir,
                          init_params = model_init_vals,
                          cohort_sub_dir = FALSE,
                          model_status = opt$model_status,
                          data_filt_list = data_filt,
                          min_iter = opt$min_iter,
                          max_iter = opt$max_iter,
                          iter_increment = opt$iter_increment,
                          diag_thresholds = diag_thresholds,
                          enable_adaptive_iter = !opt$disable_adaptive_iter)

if (!opt$dry_run) {
  cat("Hierarchical model fitted and saved successfully.\n")
  cat("Model type:", group_type, "\n")
  cat("Subjects modeled:", selected_n_subs, "\n")
} else {
  cat("Dry run completed successfully.\n")
}

cat("Script completed. Exiting explicitly.\n")
q(save = "no", status = 0)