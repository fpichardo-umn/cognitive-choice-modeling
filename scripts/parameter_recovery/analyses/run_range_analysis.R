#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(knitr)
})

# Define command line options
option_list = list(
  make_option(c("-t", "--task"), type="character", help="Task name"),
  make_option(c("-m", "--model"), type="character", help="Model name"),
  make_option(c("-g", "--group"), type="character", help="Group type"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, help="Output directory"),
  make_option(c("-r", "--render"), action="store_true", default=FALSE, help="Render RMD to HTML")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required parameters
if (is.null(opt$task) || is.null(opt$model) || is.null(opt$group)) {
  stop("Task, model, and group parameters are required.")
}

# Source helper functions
source(file.path(here::here(), "scripts", "parameter_recovery", "helper_functions_PR.R"))
source(file.path(here::here(), "scripts", "parameter_recovery", "recovery", "recovery.R"))

# Get directory structure
dirs <- setup_directories()

# Set directories
recovery_dir <- dirs$REC_SIM_DIR
ppc_dir <- dirs$PPC_DIR
output_dir <- file.path(dirs$PR_DIR, "analyses", "output")

if (!is.null(opt$outdir)) {
  output_dir <- opt$outdir
}

# Create output directory if it doesn't exist
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Determine model type
model_type <- determine_model_type(opt$model)

# Construct file paths using conventions from the existing pipeline
recovery_csv <- file.path(recovery_dir, paste0("recovery_", opt$task, "_", opt$group, "_", opt$model, ".csv"))
ppc_subject_csv <- file.path(ppc_dir, generate_bids_filename("ppc_subject_summary", opt$task, opt$group, opt$model))
ppc_model_csv <- file.path(ppc_dir, generate_bids_filename("ppc_model_summary", opt$task, opt$group, opt$model))
ppc_block_csv <- file.path(ppc_dir, generate_bids_filename("ppc_block_stats", opt$task, opt$group, opt$model))

# Define output RMD and HTML file paths
output_rmd <- file.path(output_dir, paste0("range_analysis_", opt$task, "_", opt$group, "_", opt$model, ".Rmd"))
output_html <- gsub("\\.Rmd$", ".html", output_rmd)

# Check if required files exist
if (!file.exists(recovery_csv)) {
  stop("Recovery CSV file not found: ", recovery_csv)
}

if (!file.exists(ppc_subject_csv)) {
  warning("PPC subject CSV file not found: ", ppc_subject_csv, " - Analysis will proceed without behavioral data")
  has_ppc_data <- FALSE
} else {
  has_ppc_data <- TRUE
}

# Read the template
rmd_template_file <- file.path(dirs$PR_DIR, "analyses", "range_analysis_enhanced.Rmd")
if (!file.exists(rmd_template_file)) {
  stop("RMD template file not found: ", rmd_template_file)
}

rmd_template <- readLines(rmd_template_file, warn = FALSE)
rmd_content <- paste(rmd_template, collapse = "\n")

# Initialize results file by adding parameters and file paths
modified_content <- rmd_content

# Replace template values
if (has_ppc_data) {
  modified_content <- gsub('recovery_file <- ""', paste0('recovery_file <- "', recovery_csv, '"'), modified_content)
  modified_content <- gsub('ppc_subject_file <- ""', paste0('ppc_subject_file <- "', ppc_subject_csv, '"'), modified_content)
  modified_content <- gsub('ppc_model_file <- ""', paste0('ppc_model_file <- "', ppc_model_csv, '"'), modified_content)
  
  if (file.exists(ppc_block_csv)) {
    modified_content <- gsub('ppc_block_file <- ""', paste0('ppc_block_file <- "', ppc_block_csv, '"'), modified_content)
  }
} else {
  modified_content <- gsub('recovery_file <- ""', paste0('recovery_file <- "', recovery_csv, '"'), modified_content)
}

# Add model type
modified_content <- gsub('model_type <- "RL"', paste0('model_type <- "', model_type, '"'), modified_content)

# Add task, model, and group info to title
title_replace <- paste0('title: "Enhanced Parameter Range Analysis: ', opt$task, ' - ', opt$group, ' - ', opt$model, '"')
modified_content <- sub('title: "Enhanced Parameter Range Analysis with Behavioral Integration"', title_replace, modified_content)

# Write the modified RMD
writeLines(modified_content, output_rmd)
message("Generated range analysis RMD file: ", output_rmd)

# Render to HTML if requested
if (opt$render) {
  if(requireNamespace("rmarkdown", quietly = TRUE)) {
    message("Rendering HTML report: ", output_html)
    rmarkdown::render(output_rmd, output_file = output_html)
    message("HTML report generation complete.")
  } else {
    warning("rmarkdown package not available. Cannot render HTML.")
  }
}
