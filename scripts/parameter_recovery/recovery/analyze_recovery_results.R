#!/usr/bin/env Rscript

# Script to run parameter recovery analysis with cohort and session support

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(rmarkdown)
})

# Define command line options
option_list = list(
  make_option(c("-i", "--input"), type="character", help="Input CSV file with recovery data"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output directory for RMD and HTML files"),
  make_option(c("-t", "--task"), type="character", help="Task name"),
  make_option(c("-m", "--model"), type="character", help="Model name"),
  make_option(c("-g", "--group"), type="character", help="Group type"),
  make_option(c("-r", "--render"), action="store_true", default=FALSE, help="Render to HTML"),
  make_option(c("--cohort"), type="character", required=TRUE, help="Cohort identifier"),
  make_option(c("--session"), type="character", default=NULL, help="Session identifier")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Options used:\n")
dput(opt)

# Validate required parameters
if (is.null(opt$input)) {
  stop("Input CSV file is required")
}

if (is.null(opt$task) || is.null(opt$model) || is.null(opt$group)) {
  # Try to extract information from filename
  filename_parts <- strsplit(basename(opt$input), "_")[[1]]
  
  if (length(filename_parts) >= 4) {
    if (is.null(opt$task)) opt$task <- filename_parts[2]
    if (is.null(opt$model)) opt$model <- filename_parts[3]
    if (is.null(opt$group)) opt$group <- filename_parts[4]
  } else {
    stop("Task, model, and group must be provided (could not extract from filename)")
  }
}

# Source helper functions - this loads all other helper functions
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "parameter_recovery", "recovery", "recovery.R"))

# Get directory structure
dirs <- setup_directories(opt$task)

# Set output directory
output_dir <- if (!is.null(opt$output)) {
  opt$output
} else {
  dirs$REC_SIM_DIR
}

# Make sure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define output RMD filename using BIDS-inspired format
output_rmd <- file.path(
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

# Generate RMD file
message("Generating RMD report from: ", opt$input)
message("Output RMD file: ", output_rmd)

# Use the existing function to generate the RMD file
generate_recovery_rmd(opt$input, output_rmd, opt$task, opt$model, opt$group)

# Render HTML
if (opt$render) {
  html_file <- gsub("\\.Rmd$", ".html", output_rmd)
  message("Rendering HTML report: ", html_file)
  rmarkdown::render(output_rmd, output_file = html_file)
}

message("Recovery analysis complete")
