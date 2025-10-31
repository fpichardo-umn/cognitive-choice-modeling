#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(here)
})

# Define command line options
option_list = list(
  make_option(c("--task"), type="character", help="Task name"),
  make_option(c("--model"), type="character", help="Model name"),
  make_option(c("--group"), type="character", help="Group type"),
  make_option(c("--dir"), type="character", default=NULL, help="Directory containing CSV files")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required parameters
if (is.null(opt$task) || is.null(opt$model) || is.null(opt$group)) {
  stop("Task, model, and group parameters are required.")
}

# Source helper functions
source(file.path(here::here(), "scripts", "parameter_recovery", "helper_functions_PR.R"))
source(file.path(here::here(), "scripts", "parameter_recovery", "posterior_predictive", "posterior_predictive.R"))

# Get directory structure
dirs <- setup_directories()

# Set directory
output_dir <- dirs$PPC_DIR
if (!is.null(opt$dir)) {
  output_dir <- opt$dir
}

# Determine model type based on model name
model_type <- determine_model_type(opt$model)

# Define CSV and RMD file paths using the same BIDS-like naming conventions
subject_csv <- file.path(output_dir, generate_bids_filename("ppc_subject_summary", opt$task, opt$group, opt$model))
model_csv <- file.path(output_dir, generate_bids_filename("ppc_model_summary", opt$task, opt$group, opt$model))
block_csv <- file.path(output_dir, generate_bids_filename("ppc_block_stats", opt$task, opt$group, opt$model))
rmd_file <- file.path(output_dir, paste0("ppc_analysis_", opt$task, "_", opt$group, "_", opt$model, ".Rmd"))

# Check if subject CSV file exists (required)
if (!file.exists(subject_csv)) {
  stop("Subject summary CSV file not found: ", subject_csv)
}

# Prepare input files list
input_files <- list(
  subject_csv = subject_csv,
  model_csv = model_csv,
  block_csv = block_csv
)

# Generate RMD file
message("Generating RMD report from CSV files in: ", output_dir)
message("Output RMD file: ", rmd_file)

# Use the existing function to generate the RMD file
generate_ppc_rmd(input_files, rmd_file, opt$task, opt$group, opt$model, model_type)

# Render HTML
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  html_file <- gsub("\\.Rmd$", ".html", rmd_file)
  message("Rendering HTML report: ", html_file)
  rmarkdown::render(rmd_file, output_file = html_file)
}

message("Report generation complete.")