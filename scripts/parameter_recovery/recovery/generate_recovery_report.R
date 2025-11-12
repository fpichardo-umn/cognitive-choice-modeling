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
  make_option(c("--dir"), type="character", default=NULL, help="Directory containing CSV files"),
  make_option(c("--cohort"), type="character", required=TRUE, help="Cohort identifier"),
  make_option(c("--session"), type="character", default=NULL, help="Session identifier")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Options used:\n")
dput(opt)

# Validate required parameters
if (is.null(opt$task) || is.null(opt$model) || is.null(opt$group)) {
  stop("Task, model, and group parameters are required.")
}

# Source helper functions - this loads all the other helper functions
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "parameter_recovery", "recovery", "recovery.R"))

# Get directory structure
dirs <- setup_directories(opt$task)

# Set directory
output_dir <- dirs$REC_SIM_DIR
if (!is.null(opt$dir)) {
  output_dir <- opt$dir
}

# Define CSV and RMD file paths using BIDS-inspired format
csv_file <- file.path(
  output_dir,
  generate_bids_filename(
    prefix = NULL,
    task = opt$task,
    group = opt$group,
    model = opt$model,
    cohort = opt$cohort,
    ses = opt$session,
    additional_tags = list(
      "type" = "rec",
      "desc" = "data"
    ),
    ext = "csv"
  )
)

rmd_file <- file.path(
  output_dir,
  generate_bids_filename(
    prefix = NULL,
    task = opt$task,
    group = opt$group,
    model = opt$model,
    cohort = opt$cohort,
    ses = opt$session,
    additional_tags = list(
      "type" = "rec",
      "desc" = "analysis"
    ),
    ext = "Rmd"
  )
)

# Check if CSV file exists
if (!file.exists(csv_file)) {
  stop("Recovery CSV file not found: ", csv_file)
}

# Generate RMD file
message("Generating RMD report from: ", csv_file)
message("Output RMD file: ", rmd_file)

# Use the existing function to generate the RMD file
generate_recovery_rmd(csv_file, rmd_file, opt$task, opt$model, opt$group)

# Render HTML
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  html_file <- gsub("\\.Rmd$", ".html", rmd_file)
  message("Rendering HTML report: ", html_file)
  rmarkdown::render(rmd_file, output_file = html_file)
}

message("Report generation complete.")
