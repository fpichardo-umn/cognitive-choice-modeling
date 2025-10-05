#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(parallel)
  library(dplyr)
})

# Parse command line arguments
option_list = list(
  make_option(c("-a", "--subjects"), type="character", default=NULL, 
              help="Subject indices (e.g., '40-45' or '40,42,45' or '40-45,50,52-54')"),
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="Model name"),
  make_option(c("-k", "--task"), type="character", default=NULL, 
              help="Task name"),
  make_option(c("-s", "--source"), type="character", default=NULL, 
              help="Data source"),
  make_option(c("--ses"), type="character", default=NULL, 
              help="Session identifier (optional)"),
  make_option(c("--model_status"), type="character", default=NULL,
              help="Model status (canonical/experimental/working). Default: auto-detect"),
  make_option(c("-f", "--fit_config"), type="character", default="sing", 
              help="Fit parameters config name (default: sing)"),
  make_option(c("-d", "--data_config"), type="character", default="default", 
              help="Data parameters config name (default: default)"),
  make_option(c("-t", "--type"), type="character", default="fit", 
              help="Type of stan code to run (fit, postpc, prepc) (default: fit)"),
  make_option(c("-l", "--subs_file"), type="character", default="subject_ids_complete_valid.txt", 
              help="Subs list file [Data/txt/subs/] (default: subject_ids_complete_valid.txt)"),
  make_option(c("-n", "--dry_run"), action="store_true", default=FALSE, 
              help="Perform a dry run"),
  make_option(c("-p", "--parallel"), action="store_true", default=FALSE, 
              help="Run subjects in parallel (uses multiple cores)"),
  make_option(c("-c", "--cores"), type="integer", default=4, 
              help="Number of cores to use for parallel processing (default: 4)"),
  make_option(c("--n_trials"), type = "integer", default = 120, help = "Number of trials"),
  make_option(c("--RTbound_min_ms"), type = "integer", default = 50, help = "Min RT bound in ms"),
  make_option(c("--RTbound_max_ms"), type = "integer", default = 2500, help = "Max RT bound in ms"),
  make_option(c("--rt_method"), type = "character", default = "raw", help = "RT preprocessing method"),
  make_option(c("--n_warmup"), type = "integer", default = 3000, help = "Number of warmup iterations"),
  make_option(c("--n_iter"), type = "integer", default = 15000, help = "Total number of iterations"),
  make_option(c("--n_chains"), type = "integer", default = 4, help = "Number of chains"),
  make_option(c("--adapt_delta"), type = "double", default = 0.95, help = "Stan adapt delta"),
  make_option(c("--max_treedepth"), type = "integer", default = 12, help = "Stan max tree depth"),
  make_option(c("--check_iter"), type = "integer", default = 1000, help = "Iteration interval to check convergence"),
  make_option(c("--seed"), type = "integer", help = "Random seed")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to expand array specification (e.g., "40-45,50,52-54" -> c(40,41,42,43,44,45,50,52,53,54))
expand_array_spec <- function(spec) {
  result <- integer(0)
  
  # Split by comma
  ranges <- strsplit(spec, ",")[[1]]
  
  for (range in ranges) {
    if (grepl("-", range)) {
      # It's a range
      range_parts <- as.integer(strsplit(range, "-")[[1]])
      start <- range_parts[1]
      end <- range_parts[2]
      result <- c(result, start:end)
    } else {
      # It's a single number
      result <- c(result, as.integer(range))
    }
  }
  
  return(result)
}

# Load helper functions for directory structure
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))

# Set up directories
PROJ_DIR <- get_proj_dir()
SCRIPT_DIR <- file.path(PROJ_DIR, "scripts")
SAFE_DIR <- get_safe_data_dir()

# Get task-specific directories
if (!is.null(opt$task)) {
  DATA_DIR <- get_data_dir(opt$task)
  TXT_DIR <- get_txt_dir(opt$task)
  
  # Get subjects directory with source
  SUBS_DIR <- file.path(SAFE_DIR, opt$source)
  if (!is.null(opt$ses)) {
    SUBS_DIR <- file.path(SUBS_DIR, paste0("ses-", opt$ses))
  }
  SUBS_LIST_FILE <- file.path(SUBS_DIR, opt$subs_file)
} else {
  stop("Task name is required using the -k option.")
}

# Create task-specific log directory
log_dir <- file.path(PROJ_DIR, "log_files", opt$task, "fit")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# Create log file with BIDS-style naming
# Create BIDS-style log filename
log_filename <- generate_bids_filename(
  prefix = "batch_processing",
  cohort = opt$source,
  ses = opt$ses,
  task = opt$task,
  group = "sing",
  model = opt$model,
  additional_tags = list("date" = format(Sys.time(), "%Y%m%d-%H%M%S")),
  ext = "log"
)
log_file <- file.path(log_dir, log_filename)
cat(sprintf("Log file: %s\n", log_file))

# Function to log messages with timestamps
log_message <- function(message) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  log_entry <- sprintf("%s %s", timestamp, message)
  cat(log_entry, "\n")
  cat(log_entry, "\n", file = log_file, append = TRUE)
}

# Check required arguments
if (is.null(opt$subjects) || is.null(opt$model) || is.null(opt$task) || is.null(opt$source)) {
  stop("Error: Subject indices, model name, task, and source are required.")
}

# Check subjects list file
if (!file.exists(SUBS_LIST_FILE)) {
  stop(sprintf("Subjects list file not found: %s", SUBS_LIST_FILE))
}
subject_ids <- readLines(SUBS_LIST_FILE)

# Expand subject indices
subject_indices <- expand_array_spec(opt$subjects)
total_subjects <- length(subject_indices)

# Log start of processing
log_message(sprintf("Starting batch processing for %d subjects", total_subjects))
log_message(sprintf("Model: %s, Task: %s, Source: %s", opt$model, opt$task, opt$source))
if (!is.null(opt$ses)) {
  log_message(sprintf("Session: %s", opt$ses))
}
log_message(sprintf("Fit config: %s, Data config: %s", opt$fit_config, opt$data_config))

# Validate subject indices
invalid_indices <- subject_indices[subject_indices > length(subject_ids) | subject_indices < 1]
if (length(invalid_indices) > 0) {
  log_message(sprintf("WARNING: Invalid subject indices found: %s", 
                      paste(invalid_indices, collapse = ", ")))
  subject_indices <- subject_indices[subject_indices <= length(subject_ids) & subject_indices >= 1]
  log_message(sprintf("Continuing with valid indices only (%d subjects)", length(subject_indices)))
}

# Source the config files (similar to the bash script)
fit_config_file <- file.path(SCRIPT_DIR, "configs", paste0("fit_params_", opt$fit_config, ".conf"))
data_config_file <- file.path(SCRIPT_DIR, "configs", paste0("data_params_", opt$data_config, ".conf"))

# Parse the conf files to extract parameters
# This is a simplified version - in a real implementation, you might need to adapt 
# this to match the exact structure of your conf files
parse_conf_file <- function(file_path) {
  if (!file.exists(file_path)) {
    log_message(sprintf("WARNING: Config file not found: %s", file_path))
    return(list())
  }
  
  lines <- readLines(file_path)
  params <- list()
  
  for (line in lines) {
    # Skip comments and empty lines
    if (grepl("^#", line) || grepl("^\\s*$", line)) {
      next
    }
    
    # Parse parameter assignments
    if (grepl("=", line)) {
      parts <- strsplit(line, "=")[[1]]
      param_name <- trimws(parts[1])
      param_value <- trimws(parts[2])
      
      # Remove quotes if present
      param_value <- gsub("^\"(.*)\"$", "\\1", param_value)
      param_value <- gsub("^'(.*)'$", "\\1", param_value)
      
      # Convert numeric values
      if (grepl("^[0-9]+$", param_value)) {
        param_value <- as.integer(param_value)
      } else if (grepl("^[0-9]*\\.[0-9]+$", param_value)) {
        param_value <- as.numeric(param_value)
      }
      
      params[[param_name]] <- param_value
    }
  }
  
  return(params)
}

# Parse config files
fit_params <- parse_conf_file(fit_config_file)
data_params <- parse_conf_file(data_config_file)

# Function to process a single subject
fit_single_sub <- function(index, subject_indices, subject_ids, opt, fit_params, data_params) {
  subject_id <- subject_ids[index]
  
  if (is.na(subject_id) || subject_id == "") {
    log_message(sprintf("ERROR: Could not find subject ID for index %d", index))
    return(NULL)
  }
  
  log_message(sprintf("Processing subject %s (index: %d)", subject_id, index))
  
  # Determine RT handling method (following logic from batch script)
  rt_method <- if (is.null(fit_params$rt_method)) {
    if (grepl("ddm", opt$model)) {
      "remove"
    } else {
      "force"
    }
  } else {
    fit_params$rt_method
  }
  
  # Build command arguments for Rscript
  cmd_args <- c(
    file.path(SCRIPT_DIR, "fit", "fit_single_model.R"),
    "-m", opt$model,
    "-t", opt$type,
    "-k", opt$task,
    "-s", opt$source,
    "--subid", subject_id,
    "--index", index,
    "--n_trials", as.character(opt$n_trials),
    "--RTbound_min_ms", as.character(opt$RTbound_min_ms),
    "--RTbound_max_ms", as.character(opt$RTbound_max_ms),
    "--rt_method", opt$rt_method,
    "--n_warmup", as.character(opt$n_warmup),
    "--n_iter", as.character(opt$n_iter),
    "--n_chains", as.character(opt$n_chains),
    "--adapt_delta", as.character(opt$adapt_delta),
    "--max_treedepth", as.character(opt$max_treedepth),
    "--check_iter", as.character(opt$check_iter),
    "--seed", as.character(opt$seed)
  )
  
  # Add session parameter if provided
  if (!is.null(opt$ses)) {
    cmd_args <- c(cmd_args, "--ses", opt$ses)
  }
  
  # Add model_status if provided
  if (!is.null(opt$model_status)) {
    cmd_args <- c(cmd_args, "--model_status", opt$model_status)
  }
  
  # Add --init flag if needed
  if (is.null(fit_params$init) || fit_params$init) {
    cmd_args <- c(cmd_args, "--init")
  }
  
  # Add --dry_run flag if needed
  if (opt$dry_run) {
    cmd_args <- c(cmd_args, "--dry_run")
  }
  
  # Execute the command
  if (opt$dry_run) {
    log_message(sprintf("Dry run: would execute: Rscript %s", paste(cmd_args, collapse = " ")))
    result <- 0
  } else {
    log_message(sprintf("Executing: Rscript %s", paste(cmd_args, collapse = " ")))
    result <- system2("Rscript", cmd_args)
  }
  
  if (result == 0) {
    log_message(sprintf("Successfully processed subject %s", subject_id))
    return(TRUE)
  } else {
    log_message(sprintf("ERROR: Failed to process subject %s (index: %d)", subject_id, index))
    return(FALSE)
  }
}

# Process subjects (either sequentially or in parallel)
start_time <- Sys.time()

if (opt$parallel && !opt$dry_run) {
  log_message(sprintf("Processing subjects in parallel with %d cores", opt$cores))
  
  # Create cluster
  cl <- makeCluster(min(opt$cores, total_subjects))
  
  # Export necessary variables and functions to the cluster
  clusterExport(cl, c("log_message", "SCRIPT_DIR", "opt", "fit_params", "data_params"))
  
  # Process subjects in parallel
  results <- parLapply(cl, subject_indices, function(index) {
    fit_single_sub(index, subject_indices, subject_ids, opt, fit_params, data_params)
  })
  
  # Close cluster
  stopCluster(cl)
  
  # Count successful/failed subjects
  successful <- sum(unlist(results))
  failed <- total_subjects - successful
} else {
  # Process subjects sequentially
  successful <- 0
  failed <- 0
  
  for (i in seq_along(subject_indices)) {
    index <- subject_indices[i]
    log_message(sprintf("Processing subject %d of %d", i, total_subjects))
    
    result <- fit_single_sub(index, subject_indices, subject_ids, opt, fit_params, data_params)
    
    if (!is.null(result) && result) {
      successful <- successful + 1
    } else {
      failed <- failed + 1
    }
  }
}

# Combine batch results
if (!opt$dry_run && successful > 0) {
  log_message("Combining batch of data...")
  
  # Build command arguments for combine script with proper source and session parameters
  combine_args <- c(
    file.path(SCRIPT_DIR, "fit", "combine_batch_fits.R"),
    "-k", opt$task,
    "-m", opt$model,
    "-g", "sing",
    "-s", opt$source,
    "-t", opt$type,
    "-d"
  )
  
  # Add session parameter if provided
  if (!is.null(opt$ses)) {
    combine_args <- c(combine_args, "--ses", opt$ses)
  }
  
  # Execute the command
  log_message(sprintf("Executing: Rscript %s", paste(combine_args, collapse = " ")))
  combine_result <- system2("Rscript", combine_args)
  
  if (combine_result == 0) {
    log_message("Successfully combined batch data")
  } else {
    log_message("ERROR: Failed to combine batch data")
  }
}

# Log completion
end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")
log_message(sprintf("Batch processing complete. Processed %d subjects (%d successful, %d failed) in %.2f minutes.",
                    total_subjects, successful, failed, as.numeric(elapsed)))

cat("Script completed. Exiting explicitly.\n")
