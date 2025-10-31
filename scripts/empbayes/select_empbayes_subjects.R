#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(dplyr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-k", "--task"), type="character", default=NULL, 
              help="Task name"),
  make_option(c("-s", "--source"), type="character", default=NULL, 
              help="Data source (cohort)"),
  make_option(c("--ses"), type="character", default=NULL, 
              help="Session identifier (optional)"),
  make_option(c("--n_hier"), type="integer", default=NULL, 
              help="Number of subjects for hierarchical model"),
  make_option(c("--n_trials"), type="integer", default=120, 
              help="Minimum number of trials required"),
  make_option(c("--RTbound_min_ms"), type="integer", default=50, 
              help="RT minimum bound in milliseconds"),
  make_option(c("--RTbound_max_ms"), type="integer", default=4000, 
              help="RT maximum bound in milliseconds"),
  make_option(c("--rt_method"), type="character", default="mark", 
              help="RT preprocessing method"),
  make_option(c("--seed"), type="integer", default=29518, 
              help="Random seed for sampling"),
  make_option(c("-l", "--subs_file"), type="character", default="subject_ids_complete_valid.txt", 
              help="Subs list file [Data/raw/COHORT/ses-SES/] (default: subject_ids_complete_valid.txt)"),
  make_option(c("--hier_subs_file"), type="character", default=NULL, 
              help="Pre-specified subject list file (skips random sampling)"),
  make_option(c("--dry_run"), action="store_true", default=FALSE, 
              help="Perform a dry run")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$task) || is.null(opt$source)) {
  stop("Task and source are required. Use -k and -s options.")
}

if (is.null(opt$n_hier) && is.null(opt$hier_subs_file)) {
  stop("Either --n_hier or --hier_subs_file must be provided.")
}

# Load helper functions
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))

# Set random seed
set.seed(opt$seed)

# Set up output directory
empbayes_subs_dir <- file.path(get_proj_dir(), "Outputs", opt$task, "empbayes", "subs")
ensure_dir_exists(empbayes_subs_dir)

# Generate output filename using BIDS naming
output_filename <- generate_bids_filename(
  prefix = NULL,
  task = opt$task,
  group = "emp",
  cohort = opt$source,
  ses = opt$ses,
  model = "gen",
  additional_tags = list("desc" = "hier_subs"),
  ext = "txt"
)
output_file <- file.path(empbayes_subs_dir, output_filename)

# Get task-specific directories
if (!is.null(opt$task)) {
  SAFE_DIR <- get_safe_data_dir()
  # Get subjects directory with source
  SUBS_DIR <- file.path(SAFE_DIR, opt$source)
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

cat("Empirical Bayes Subject Selection\n")
cat("==================================\n\n")

# If pre-specified file provided, validate and copy it
if (!is.null(opt$hier_subs_file)) {
  if (!file.exists(opt$hier_subs_file)) {
    stop("Specified hier_subs_file does not exist: ", opt$hier_subs_file)
  }
  
  cat("Using pre-specified subject list from:", opt$hier_subs_file, "\n")
  
  if (!opt$dry_run) {
    file.copy(opt$hier_subs_file, output_file, overwrite = TRUE)
    selected_subjects <- readLines(output_file)
    cat("Copied", length(selected_subjects), "subjects to:", output_file, "\n")
  } else {
    cat("Dry run: Would copy file to:", output_file, "\n")
  }
  
  quit(status = 0)
}

# Otherwise, perform random sampling
cat("Loading data...\n")

if (!opt$dry_run) {
  # Load data
  all_data <- load_data(opt$task, opt$source, opt$ses)
  
  # Filter based on main subs file
  all_data = all_data[all_data$subjID %in% subject_ids,]
  
  # Get unique subjects
  all_subjects <- unique(all_data$subjID)
  cat("Total subjects in dataset:", length(all_subjects), "\n")
  
  # Apply data quality filters to identify eligible subjects
  cat("\nApplying data quality filters...\n")
  cat("  Minimum trials:", opt$n_trials, "\n")
  cat("  RT bounds:", opt$RTbound_min_ms, "-", opt$RTbound_max_ms, "ms\n")
  cat("  RT method:", opt$rt_method, "\n\n")
  
  eligible_subjects <- c()
  
  for (subj in all_subjects) {
    subj_data <- all_data[all_data$subjID == subj, ]
    
    # Check minimum trial count
    if (nrow(subj_data) < opt$n_trials) {
      next
    }
    
    # Apply RT filtering based on method
    if (opt$rt_method == "remove") {
      valid_trials <- subj_data$RT >= (opt$RTbound_min_ms / 1000) & 
                      subj_data$RT <= (opt$RTbound_max_ms / 1000)
      remaining_trials <- sum(valid_trials)
      
      if (remaining_trials < opt$n_trials) {
        next
      }
    }
    
    eligible_subjects <- c(eligible_subjects, subj)
  }
  
  cat("Eligible subjects after filtering:", length(eligible_subjects), "\n")
  
  if (length(eligible_subjects) < opt$n_hier) {
    stop("Not enough eligible subjects. Requested: ", opt$n_hier, 
         ", Available: ", length(eligible_subjects))
  }
  
  # Random sample
  selected_subjects <- sample(eligible_subjects, opt$n_hier, replace = FALSE)
  
  cat("\nRandomly selected", opt$n_hier, "subjects (seed:", opt$seed, ")\n")
  
  # Save to file
  writeLines(as.character(selected_subjects), output_file)
  cat("Saved subject list to:", output_file, "\n")
  
  # Print first few subjects
  cat("\nFirst 10 selected subjects:\n")
  print(head(selected_subjects, 10))
  
} else {
  cat("Dry run: Would load data from load_data(", opt$task, ", ", opt$source, ", ", opt$ses, ")\n")
  cat("Dry run: Would apply filters:\n")
  cat("  n_trials >=", opt$n_trials, "\n")
  cat("  RT bounds:", opt$RTbound_min_ms, "-", opt$RTbound_max_ms, "ms\n")
  cat("  RT method:", opt$rt_method, "\n")
  cat("Dry run: Would randomly sample", opt$n_hier, "subjects\n")
  cat("Dry run: Would save to:", output_file, "\n")
}

cat("\nSubject selection complete.\n")
