#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(cmdstanr)
  library(posterior)
  library(dplyr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="Model name"),
  make_option(c("-k", "--task"), type="character", default=NULL, 
              help="Task name"),
  make_option(c("-s", "--source"), type="character", default=NULL, 
              help="Data source (cohort)"),
  make_option(c("--ses"), type="character", default=NULL, 
              help="Session identifier (optional)"),
  make_option(c("--model_status"), type="character", default=NULL,
              help="Model status (canonical/experimental/working). Default: auto-detect"),
  make_option(c("-d", "--data"), type="character", default=NULL, 
              help="Comma-separated list of data to extract"),
  make_option(c("-p", "--params"), type="character", default=NULL, 
              help="Comma-separated list of model parameters"),
  make_option(c("--n_trials"), type="integer", default=120, 
              help="Number of trials"),
  make_option(c("--RTbound_min_ms"), type="integer", default=50, 
              help="RT minimum bound in milliseconds"),
  make_option(c("--RTbound_max_ms"), type="integer", default=2500, 
              help="RT maximum bound in milliseconds"),
  make_option(c("--rt_method"), type="character", default="mark", 
              help="RT method"),
  make_option(c("--n_warmup"), type="integer", default=3000, 
              help="Number of warmup iterations"),
  make_option(c("--n_iter"), type="integer", default=15000, 
              help="Number of iterations"),
  make_option(c("--n_chains"), type="integer", default=4, 
              help="Number of chains"),
  make_option(c("--adapt_delta"), type="double", default=0.95, 
              help="Adapt delta"),
  make_option(c("--max_treedepth"), type="integer", default=12, 
              help="Max tree depth"),
  make_option(c("--seed"), type="integer", default=29518, 
              help="Set seed"),
  make_option(c("--min_valid_rt_pct"), type="double", default=0.7, help="Minimum percent valid RT"),
  make_option(c("--check_iter"), type="integer", default=5000, 
              help="Iteration interval for checkpoint runs"),
  make_option(c("--hier_subs_file"), type="character", default=NULL, 
              help="Path to hierarchical subject list file"),
  make_option(c("--fitting_method"), type="character", default="mcmc",
              help="Fitting method: mcmc or pathfinder (default: mcmc)"),
  make_option(c("--pf_num_paths"), type="integer", default=4,
              help="Pathfinder: number of paths (default: 4)"),
  make_option(c("--pf_draws"), type="integer", default=1000,
              help="Pathfinder: number of final draws after PSIS (default: 1000)"),
  make_option(c("--pf_single_path_draws"), type="integer", default=250,
              help="Pathfinder: draws per single path (default: 250)"),
  make_option(c("--dry_run"), action="store_true", default=FALSE, 
              help="Perform a dry run")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$model) || is.null(opt$task) || is.null(opt$source)) {
  stop("Model, task, and source are required. Use -m, -k, and -s options.")
}

# Load helper functions
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))

# Set random seed
set.seed(opt$seed)

# Get model defaults
model_defaults <- get_model_defaults(opt$task)

# Build model name
group_type <- "hier"
full_model_name <- paste(opt$task, group_type, opt$model, sep="_")

if (!full_model_name %in% names(model_defaults)) {
  cat("Requested model:", full_model_name, "\n")
  cat("Available models:", names(model_defaults), "\n")
  stop("Unrecognized hierarchical model. Please check the model name.")
}

# Extract data and parameter specifications
data_to_extract <- if (!is.null(opt$data)) {
  strsplit(opt$data, ",")[[1]]
} else {
  model_defaults[[full_model_name]]$data
}

model_params <- if (!is.null(opt$params)) {
  strsplit(opt$params, ",")[[1]]
} else {
  model_defaults[[full_model_name]]$params
}

cat("Empirical Bayes Hierarchical Model Fitting\n")
cat("==========================================\n\n")
cat("Model:", full_model_name, "\n")
cat("Task:", opt$task, "\n")
cat("Source:", opt$source, "\n")
if (!is.null(opt$ses)) cat("Session:", opt$ses, "\n")
cat("\n")

# Determine subject list file
if (is.null(opt$hier_subs_file)) {
  # Use default location
  empbayes_subs_dir <- file.path(get_proj_dir(), "Outputs", opt$task, "empbayes", "subs")
  hier_subs_filename <- generate_bids_filename(
    prefix = NULL,
    task = opt$task,
    group = "emp",
    cohort = opt$source,
    ses = opt$ses,
    model = "gen",
    additional_tags = list("desc" = "hier_subs"),
    ext = "txt"
  )
  opt$hier_subs_file <- file.path(empbayes_subs_dir, hier_subs_filename)
}

# Check that subject list exists
if (!file.exists(opt$hier_subs_file)) {
  stop("Hierarchical subject list not found: ", opt$hier_subs_file, 
       "\nRun select_empbayes_subjects.R first.")
} else {
  cat("Using this file for subjects:\n")
  cat(opt$hier_subs_file)
}

# Load subject list
hier_subjects <- readLines(opt$hier_subs_file)
n_hier_subs <- length(hier_subjects)

cat("Loading", n_hier_subs, "subjects from:", basename(opt$hier_subs_file), "\n\n")

# Load and prepare data
if (!opt$dry_run) {
  cat("Loading data...\n")
  all_data <- load_data(opt$task, opt$source, opt$ses)
  
  # Filter to hierarchical subjects
  hierarchical_data <- all_data[all_data$subjID %in% hier_subjects, ]
  
  # Verify we have all subjects
  found_subjects <- unique(hierarchical_data$subjID)
  if (length(found_subjects) != n_hier_subs) {
    missing <- setdiff(hier_subjects, found_subjects)
    warning("Some subjects from list not found in data. Missing: ", 
            paste(missing, collapse=", "))
  }
  
  cat("Extracted data for", length(found_subjects), "subjects\n")
  
  # Extract data for hierarchical model
  data_list <- extract_sample_data(
    hierarchical_data, 
    data_to_extract, 
    task = opt$task,
    n_subs = length(found_subjects), 
    n_trials = opt$n_trials, 
    RTbound_min_ms = opt$RTbound_min_ms, 
    RTbound_max_ms = opt$RTbound_max_ms,
    RTbound_reject_min_ms = opt$RTbound_min_ms + 20, 
    RTbound_reject_max_ms = opt$RTbound_max_ms, 
    rt_method = opt$rt_method, 
    minrt_ep_ms = 0, min_valid_rt_pct = opt$min_valid_rt_pct
  )
  
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
  
  # Handle entropy data if present
  if (sum(grepl("entropy", names(data_list))) > 0) {
    cat("Processing entropy measures...\n")
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
  cat("Dry run: Would load data using load_data(", opt$task, ", ", opt$source, ", ", opt$ses, ")\n")
  cat("Dry run: Would filter to", n_hier_subs, "subjects\n")
  cat("Dry run: Data to extract:", paste(data_to_extract, collapse=", "), "\n")
  data_list <- NULL
  found_subjects <- hier_subjects
}

# Set output directory for hierarchical models
output_dir <- get_empbayes_output_dir(opt$task, "hierarchical", opt$source)

cat("\nFitting hierarchical model...\n")
cat("Output directory:", output_dir, "\n\n")

# Fit hierarchical model
fit <- fit_and_save_model(
  task = opt$task,
  cohort = opt$source,
  ses = opt$ses,
  group_type = group_type,
  model_name = opt$model,
  model_type = "fit",
  data_list = data_list,
  n_subs = length(found_subjects),
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
  init_params = NULL,
  cohort_sub_dir = FALSE,
  model_status = opt$model_status,
  data_filt_list = data_filt,
  fitting_method = opt$fitting_method,
  pf_num_paths = opt$pf_num_paths,
  pf_draws = opt$pf_draws,
  pf_single_path_draws = opt$pf_single_path_draws
)

if (!opt$dry_run) {
  cat("\nHierarchical model fitted and saved successfully.\n")
  cat("Subjects modeled:", length(found_subjects), "\n")
  cat("\nSubject list saved in fit object for tracking.\n")
} else {
  cat("\nDry run completed successfully.\n")
}
