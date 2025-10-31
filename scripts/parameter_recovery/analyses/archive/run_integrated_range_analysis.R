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
  make_option(c("--recovery"), type="character", help="Recovery CSV file path"),
  make_option(c("--ppc_subject"), type="character", help="PPC subject summary CSV file path"),
  make_option(c("--ppc_model"), type="character", help="PPC model summary CSV file path"),
  make_option(c("--ppc_block"), type="character", help="PPC block stats CSV file path (optional)", default=NULL),
  make_option(c("--task"), type="character", help="Task name"),
  make_option(c("--model"), type="character", help="Model name"),
  make_option(c("--group"), type="character", help="Group type"),
  make_option(c("--outdir"), type="character", default=NULL, help="Output directory"),
  make_option(c("--render"), action="store_true", default=FALSE, help="Render RMD to HTML")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required parameters
if (is.null(opt$recovery) || is.null(opt$ppc_subject) || is.null(opt$ppc_model) || 
    is.null(opt$task) || is.null(opt$model) || is.null(opt$group)) {
  stop("Recovery file, PPC files, task, model, and group parameters are required.")
}

# Source helper functions
source(file.path(here::here(), "scripts", "parameter_recovery", "helper_functions_PR.R"))
source(file.path(here::here(), "scripts", "parameter_recovery", "recovery", "recovery.R"))

# Get directory structure
dirs <- setup_directories()

# Set output directory
output_dir <- file.path(dirs$PR_DIR, "analyses", "results")
if (!is.null(opt$outdir)) {
  output_dir <- opt$outdir
}

# Make sure output directory exists
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Determine model type
determine_model_type <- function(model_name) {
  # Determine model type based on model name
  if(grepl("_ddm|_ssm", model_name, ignore.case = TRUE)) {
    return("RL_SSM")  # Combined model type
  } else if(grepl("^ddm$|^ssm$", model_name, ignore.case = TRUE)) {
    return("SSM")  # Sequential Sampling Model
  } else {
    return("RL")  # Reinforcement Learning model
  }
}

model_type <- determine_model_type(opt$model)

# Define output RMD file path
output_rmd <- file.path(output_dir, paste0("range_analysis_integrated_", opt$task, "_", opt$model, "_", opt$group, ".Rmd"))

# Check if input files exist
for (file_path in c(opt$recovery, opt$ppc_subject, opt$ppc_model)) {
  if (!file.exists(file_path)) {
    stop("Input file does not exist: ", file_path)
  }
}

# Generate the RMD file with appropriate placeholders replaced
rmd_template_path <- file.path(here::here(), "scripts", "parameter_recovery", "analyses", "range_analysis_enhanced_template.Rmd")
if(!file.exists(rmd_template_path)) {
  stop("RMD template not found at: ", rmd_template_path)
}

# Read the template
template_content <- readLines(rmd_template_path, warn = FALSE)
template_text <- paste(template_content, collapse = "\n")

# Replace placeholders
rmd_content <- gsub("\\{\\{TASK\\}\\}", opt$task, template_text)
rmd_content <- gsub("\\{\\{MODEL\\}\\}", opt$model, rmd_content)
rmd_content <- gsub("\\{\\{GROUP\\}\\}", opt$group, rmd_content)
rmd_content <- gsub("\\{\\{MODEL_TYPE\\}\\}", model_type, rmd_content)
rmd_content <- gsub("\\{\\{RECOVERY_FILE\\}\\}", opt$recovery, rmd_content)
rmd_content <- gsub("\\{\\{PPC_SUBJECT_FILE\\}\\}", opt$ppc_subject, rmd_content)
rmd_content <- gsub("\\{\\{PPC_MODEL_FILE\\}\\}", opt$ppc_model, rmd_content)
rmd_content <- gsub("\\{\\{PPC_BLOCK_FILE\\}\\}", ifelse(is.null(opt$ppc_block), "", opt$ppc_block), rmd_content)

# Write to output file
writeLines(rmd_content, output_rmd)
message("Integrated range analysis RMD file generated: ", output_rmd)

# Render HTML if requested
if (opt$render && requireNamespace("rmarkdown", quietly = TRUE)) {
  html_file <- gsub("\\.Rmd$", ".html", output_rmd)
  message("Rendering RMD to HTML...")
  rmarkdown::render(output_rmd, output_file = html_file)
  message("HTML file generated: ", html_file)
}

message("\nAnalysis completed. To view the results, open the generated file in R Studio or view the HTML if rendered.")
