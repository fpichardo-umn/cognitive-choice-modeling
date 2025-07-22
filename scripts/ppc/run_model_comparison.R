#!/usr/bin/env Rscript

#' Run Model Comparison
#' @description Compare multiple models based on PPC and LOO metrics

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Define command line options
option_list = list(
  make_option(c("-m", "--models"), type="character", 
              help="Comma-separated list of model names to compare"),
  make_option(c("-k", "--task"), type="character", default="igt_mod", 
              help="Task name"),
  make_option(c("-g", "--group"), type="character", default="batch_001", 
              help="Group name"),
  make_option(c("-s", "--stats_dir"), type="character", default=NULL, 
              help="Directory containing PPC statistics (optional)"),
  make_option(c("-l", "--loglik_dir"), type="character", default=NULL, 
              help="Directory containing log-likelihood data (optional)"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory (optional)"),
  make_option(c("-f", "--force"), action="store_true", default=FALSE, 
              help="Force overwrite of existing report")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$models)) {
  stop("List of models is required. Use --models flag.")
}

# Split models string into a vector of model names
models_to_compare <- strsplit(opt$models, ",")[[1]]
if (length(models_to_compare) < 2) {
  stop("At least two models are required for comparison.")
}

# Source helper functions
script_dir <- file.path(here::here(), "scripts")
source(file.path(script_dir, "parameter_recovery", "helper_functions_PR.R"))
source(file.path(script_dir, "helpers", "helper_functions_cmdSR.R"))
source(file.path(script_dir, "ppc", "model_comparison_functions.R"))

# Get directories
dirs <- setup_directories(create_missing = TRUE)
PPC_DIR <- file.path(dirs$DATA_DIR, "ppc")
PPC_STATS_DIR <- file.path(PPC_DIR, "stats")
PPC_LOGLIK_DIR <- file.path(PPC_DIR, "loglik")
PPC_REPORTS_DIR <- file.path(PPC_DIR, "reports")

# Create directories if they don't exist
for (dir in c(PPC_DIR, PPC_STATS_DIR, PPC_LOGLIK_DIR, PPC_REPORTS_DIR)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# Set directories
stats_dir <- if(!is.null(opt$stats_dir)) opt$stats_dir else PPC_STATS_DIR
loglik_dir <- if(!is.null(opt$loglik_dir)) opt$loglik_dir else PPC_LOGLIK_DIR
output_dir <- if(!is.null(opt$output_dir)) opt$output_dir else PPC_REPORTS_DIR

# Generate output filename
output_filename <- file.path(output_dir, 
                            paste0("model_comparison_", opt$task, "_", opt$group, ".html"))

# Check if output already exists
if (file.exists(output_filename) && !opt$force) {
  stop("Output file already exists. Use --force to overwrite: ", output_filename)
}

# Load PPC statistics for all models
message("Loading PPC statistics for models: ", paste(models_to_compare, collapse = ", "))
model_stats_files <- list()

for (model in models_to_compare) {
  # Construct expected filename
  stats_file <- file.path(stats_dir, 
                        generate_bids_filename(
                          prefix = "ppc_summary",
                          task = opt$task,
                          group = opt$group,
                          model = model,
                          ext = "csv"
                        ))
  
  if (file.exists(stats_file)) {
    model_stats_files[[model]] <- stats_file
    message("Found statistics file for model: ", model)
  } else {
    warning("Statistics file not found for model: ", model)
  }
}

# Check if we have enough models
if (length(model_stats_files) < 2) {
  stop("Not enough models with statistics found. Needed at least 2, found ", 
       length(model_stats_files))
}

# Load PPC statistics
ppc_stats_list <- list()
for (model in names(model_stats_files)) {
  file_path <- model_stats_files[[model]]
  ppc_stats_list[[model]] <- read.csv(file_path)
  ppc_stats_list[[model]]$model <- model  # Add model name column
  ppc_stats_list[[model]]$filename <- file_path  # Add file path
}

# Load log-likelihood files if available
loglik_list <- list()
for (model in models_to_compare) {
  # Construct expected filename
  loglik_file <- file.path(loglik_dir, 
                         generate_bids_filename(
                           prefix = "ppc_loglik",
                           task = opt$task,
                           group = opt$group,
                           model = model,
                           ext = "rds"
                         ))
  
  if (file.exists(loglik_file)) {
    loglik_list[[model]] <- readRDS(loglik_file)
    message("Found log-likelihood file for model: ", model)
  }
}

# Compare models using PPC statistics
message("Comparing models based on PPC performance...")
ppc_comparison <- compare_models_ppc(ppc_stats_list)

# Generate comparison plots
message("Generating comparison plots...")
comparison_plots <- list()
if (nrow(ppc_comparison) > 1) {
  comparison_plots$extreme_rate <- generate_comparison_plots(ppc_comparison, "extreme_rate")
  comparison_plots$ppp_dist <- generate_comparison_plots(ppc_comparison, "ppp_dist")
  comparison_plots$component <- generate_comparison_plots(ppc_comparison, "component")
}

# Compare models using LOO if available
loo_comparison <- NULL
if (length(loglik_list) >= 2) {
  message("Comparing models based on LOO-CV...")
  loo_comparison <- compare_models_loo(loglik_list)
}

# Combine results
comparison_results <- list(
  model_ppc_list = ppc_stats_list,
  ppc_comparison = ppc_comparison,
  loo_comparison = loo_comparison,
  plots = comparison_plots
)

# Generate comparison report
message("Generating comparison report...")
report_file <- generate_comparison_report(comparison_results, output_filename)

message("Model comparison complete. Report saved to: ", report_file)

# Save comparison results 
comparison_rds <- file.path(output_dir, 
                          paste0("model_comparison_", opt$task, "_", opt$group, ".rds"))
saveRDS(comparison_results, comparison_rds)
message("Comparison results saved to: ", comparison_rds)
