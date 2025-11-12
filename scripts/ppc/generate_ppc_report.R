#!/usr/bin/env Rscript

#' Generate Report for Posterior Predictive Checks
#' @description Create HTML report with PPC statistics and visualizations

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(dplyr)
  library(ggplot2)
  library(knitr)
  library(rmarkdown)
})

# Define command line options
option_list = list(
  make_option(c("-m", "--model"), type="character", default=NULL, help="Model name (e.g., ev, pvldelta)"),
  make_option(c("-k", "--task"), type="character", default=NULL, help="Task name (igt, igt_mod)"),
  make_option(c("-c", "--cohort"), type="character", default=NULL, help="Cohort identifier"),
  make_option(c("-g", "--group"), type="character", default="batch_001", help="Group identifier for file naming"),
  make_option(c("--ses"), type="character", default=NULL, help="Session identifier"),
  make_option(c("-s", "--stats_file"), type="character", default=NULL, 
              help="Path to statistics file (optional, built from task, group, model if not provided)"),
  make_option(c("-l", "--loglik_file"), type="character", default=NULL, 
              help="Path to log-likelihood file (optional)"),
  make_option(c("-e", "--exclude_file"), type="character", default=NULL, 
              help="Path to file with subject IDs to exclude (optional)"),
  make_option(c("-t", "--template"), type="character", default=NULL,
              help="Path to custom Rmd template (optional)"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory (optional)"),
  make_option(c("-f", "--force"), action="store_true", default=FALSE,
              help="Force regeneration of report even if it exists"),
  make_option(c("-r", "--render"), action="store_true", default=FALSE,
              help="Render the Rmd file to HTML after generating it")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$model) || is.null(opt$task) || is.null(opt$cohort)) {
  stop("Model name, task name, and cohort are all required.")
}

# Import helper modules - use unified helper entry point
script_dir <- file.path(here::here(), "scripts")
source(file.path(script_dir, "ppc", "helpers", "helper_functions_ppc.R"))
source(file.path(here::here(), "scripts", "ppc", "helpers", "helper_ppc_report.R"))

# Ensure PPC directories exist
directory_info <- ensure_ppc_dirs(opt$task, opt$cohort, opt$ses)

# Set output directory
output_dir <- if(!is.null(opt$output_dir)) opt$output_dir else get_ppc_reports_dir(opt$task, opt$cohort, opt$ses)

# Build statistics file path if not provided
if (is.null(opt$stats_file)) {
  # Use helper function to get the PPC statistics summary file path
  stats_file <- file.path(get_ppc_stats_dir(opt$task, opt$cohort, opt$ses),
                         generate_bids_filename("ppc_summary", opt$task, opt$group, opt$model, 
                                               "csv", cohort = opt$cohort, ses = opt$ses))
  message("Using statistics file: ", stats_file)
} else {
  stats_file <- opt$stats_file
}

# Load statistics file with error handling
stats_data <- load_csv_safe(stats_file, "PPC statistics")
if (is.null(stats_data)) {
  stop("Failed to load statistics from: ", stats_file)
}

# Build log-likelihood file path if specified
loglik_file <- NULL
if (!is.null(opt$loglik_file)) {
  loglik_file <- opt$loglik_file
} else {
  # Try to find log-likelihood file with standard path
  potential_loglik_file <- get_ppc_loglik_file_path(opt$task, opt$model, opt$group, opt$cohort, opt$ses)
  if (file.exists(potential_loglik_file)) {
    loglik_file <- potential_loglik_file
    message("Found log-likelihood file: ", loglik_file)
  }
}

# Generate the report
message("Generating PPC report...")
report_files <- generate_ppc_report(
  task_name = opt$task,
  model_name = opt$model,
  group_name = opt$group,
  cohort = opt$cohort,
  session = opt$ses,
  stats_file = stats_file,
  loglik_file = loglik_file,
  template_file = opt$template,
  exclude_file = opt$exclude_file,
  force = opt$force,
  render = opt$render
)

message("Report generation complete.")
message("Rmd file saved to: ", report_files$rmd_file)
if (!is.null(report_files$html_file)) {
  message("HTML file saved to: ", report_files$html_file)
}

