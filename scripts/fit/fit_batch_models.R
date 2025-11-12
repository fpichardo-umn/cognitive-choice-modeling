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
  make_option(c("-o", "--overwrite"), action="store_true", default=FALSE, 
              help="Overwrite existing fits (default: skip if fit already exists)"),
  make_option(c("-p", "--parallel"), action="store_true", default=FALSE, 
              help="Run subjects in parallel (uses multiple cores)"),
  make_option(c("-c", "--cores"), type="integer", default=4, 
              help="Number of cores to use for parallel processing (default: 4)"),
  make_option(c("--no_retry"), action="store_true", default=FALSE,
              help="Disable automatic retry of failed subjects"),
  make_option(c("--keep_checkpoints"), action="store_true", default=FALSE,
              help="Keep checkpoint files after successful fits (for debugging)"),
  make_option(c("--n_trials"), type = "integer", default = 120, help = "Number of trials"),
  make_option(c("--RTbound_min_ms"), type = "integer", default = 50, help = "Min RT bound in ms"),
  make_option(c("--RTbound_max_ms"), type = "integer", default = 4000, help = "Max RT bound in ms"),
  make_option(c("--rt_method"), type = "character", default = "raw", help = "RT preprocessing method"),
  make_option(c("--n_warmup"), type = "integer", default = 3000, help = "Number of warmup iterations"),
  make_option(c("--n_iter"), type = "integer", default = 15000, help = "Total number of iterations"),
  make_option(c("--n_chains"), type = "integer", default = 4, help = "Number of chains"),
  make_option(c("--adapt_delta"), type = "double", default = 0.95, help = "Stan adapt delta"),
  make_option(c("--max_treedepth"), type = "integer", default = 12, help = "Stan max tree depth"),
  make_option(c("--check_iter"), type = "integer", default = 1000, help = "Iteration interval to check convergence"),
  make_option(c("--seed"), type = "integer", help = "Random seed"),
  make_option(c("--min_iter"), type = "integer", default = NULL, 
              help = "Minimum post-warmup iterations for adaptive fitting (enables adaptive mode)"),
  make_option(c("--max_iter"), type = "integer", default = NULL,
              help = "Maximum post-warmup iterations for adaptive fitting"),
  make_option(c("--iter_increment"), type = "integer", default = 1000,
              help = "Iterations to add at each diagnostic check (default: 1000)"),
  make_option(c("--target_rhat"), type = "double", default = 1.01,
              help = "Target Rhat threshold (default: 1.01)"),
  make_option(c("--target_ess_bulk"), type = "integer", default = 400,
              help = "Target ESS bulk threshold (default: 400)"),
  make_option(c("--target_ess_tail"), type = "integer", default = 400,
              help = "Target ESS tail threshold (default: 400)"),
  make_option(c("--disable_adaptive_iter"), action = "store_true", default = FALSE,
              help = "Disable adaptive iteration feature (use fixed n_iter)")
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
log_message(sprintf("Overwrite mode: %s", if(opt$overwrite) "ON (re-fit all)" else "OFF (skip existing fits)"))
log_message(sprintf("Retry mode: %s", if(opt$no_retry) "DISABLED" else "ENABLED"))
log_message(sprintf("Checkpoint cleanup: %s", if(opt$keep_checkpoints) "DISABLED (keep all)" else "ENABLED (clean successful)"))

# Validate subject indices
invalid_indices <- subject_indices[subject_indices > length(subject_ids) | subject_indices < 1]
if (length(invalid_indices) > 0) {
  log_message(sprintf("WARNING: Invalid subject indices found: %s", 
                      paste(invalid_indices, collapse = ", ")))
  subject_indices <- subject_indices[subject_indices <= length(subject_ids) & subject_indices >= 1]
  log_message(sprintf("Continuing with valid indices only (%d subjects)", length(subject_indices)))
}

# Source the config files
fit_config_file <- file.path(SCRIPT_DIR, "configs", paste0("fit_params_", opt$fit_config, ".conf"))
data_config_file <- file.path(SCRIPT_DIR, "configs", paste0("data_params_", opt$data_config, ".conf"))

# Parse the conf files to extract parameters
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

# Function to get checkpoint file path for a subject
get_checkpoint_path <- function(subject_id, opt) {
  model_status <- if (!is.null(opt$model_status)) opt$model_status else "canonical"
  fit_dir <- file.path(get_safe_data_dir(), opt$source)
  if (!is.null(opt$ses)) {
    fit_dir <- file.path(fit_dir, paste0("ses-", opt$ses))
  }
  fit_dir <- file.path(fit_dir, "fits", opt$type, opt$task, model_status, opt$model)
  
  fit_filename <- sprintf("sub-%s_task-%s_model-%s_fit.rds", subject_id, opt$task, opt$model)
  fit_path <- file.path(fit_dir, fit_filename)
  checkpoint_path <- paste0(tools::file_path_sans_ext(fit_path), "_checkpoint.rds")
  
  return(checkpoint_path)
}

# Check for existing checkpoint files at start
check_existing_checkpoints <- function(subject_indices, subject_ids, opt) {
  existing_checkpoints <- c()
  
  for (index in subject_indices) {
    subject_id <- subject_ids[index]
    if (is.na(subject_id) || subject_id == "") next
    
    checkpoint_path <- get_checkpoint_path(subject_id, opt)
    if (file.exists(checkpoint_path)) {
      existing_checkpoints <- c(existing_checkpoints, subject_id)
    }
  }
  
  return(existing_checkpoints)
}

# Check for existing checkpoints
existing_checkpoints <- check_existing_checkpoints(subject_indices, subject_ids, opt)
if (length(existing_checkpoints) > 0) {
  log_message(sprintf("WARNING: Found %d checkpoint files from previous runs:", length(existing_checkpoints)))
  log_message(sprintf("  Subjects: %s", paste(existing_checkpoints, collapse = ", ")))
  log_message("These may indicate incomplete fits that need attention.")
}

# Function to process a single subject
fit_single_sub <- function(index, subject_indices, subject_ids, opt, fit_params, data_params) {
  subject_id <- subject_ids[index]
  
  if (is.na(subject_id) || subject_id == "") {
    log_message(sprintf("ERROR: Could not find subject ID for index %d", index))
    return(list(success = FALSE, subject_id = NA))
  }
  
  # Construct expected output path to check if fit already exists
  model_status <- if (!is.null(opt$model_status)) opt$model_status else "canonical"
  fit_dir <- file.path(get_safe_data_dir(), opt$source)
  if (!is.null(opt$ses)) {
    fit_dir <- file.path(fit_dir, paste0("ses-", opt$ses))
  }
  fit_dir <- file.path(fit_dir, "fits", opt$type, opt$task, model_status, opt$model)
  
  # Expected output filename (BIDS-style)
  fit_filename <- sprintf("sub-%s_task-%s_model-%s_fit.rds", subject_id, opt$task, opt$model)
  fit_path <- file.path(fit_dir, fit_filename)
  
  # Check if fit already exists and skip if overwrite is FALSE
  if (file.exists(fit_path) && !opt$overwrite) {
    log_message(sprintf("Skipping subject %s (index: %d) - fit already exists at: %s", 
                        subject_id, index, fit_path))
    return(list(success = TRUE, subject_id = subject_id))  # Return success since the fit exists
  }
  
  if (file.exists(fit_path) && opt$overwrite) {
    log_message(sprintf("Overwriting existing fit for subject %s (index: %d)", subject_id, index))
  } else {
    log_message(sprintf("Processing subject %s (index: %d)", subject_id, index))
  }
  
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
  
  # Add adaptive iteration parameters if provided
  if (!is.null(opt$min_iter)) {
    cmd_args <- c(cmd_args, "--min_iter", as.character(opt$min_iter))
  }
  if (!is.null(opt$max_iter)) {
    cmd_args <- c(cmd_args, "--max_iter", as.character(opt$max_iter))
  }
  if (!is.null(opt$iter_increment) && opt$iter_increment != 1000) {
    cmd_args <- c(cmd_args, "--iter_increment", as.character(opt$iter_increment))
  }
  if (!is.null(opt$target_rhat) && opt$target_rhat != 1.01) {
    cmd_args <- c(cmd_args, "--target_rhat", as.character(opt$target_rhat))
  }
  if (!is.null(opt$target_ess_bulk) && opt$target_ess_bulk != 400) {
    cmd_args <- c(cmd_args, "--target_ess_bulk", as.character(opt$target_ess_bulk))
  }
  if (!is.null(opt$target_ess_tail) && opt$target_ess_tail != 400) {
    cmd_args <- c(cmd_args, "--target_ess_tail", as.character(opt$target_ess_tail))
  }
  if (opt$disable_adaptive_iter) {
    cmd_args <- c(cmd_args, "--disable_adaptive_iter")
  }
  
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
    return(list(success = TRUE, subject_id = subject_id))
  } else {
    log_message(sprintf("ERROR: Failed to process subject %s (index: %d)", subject_id, index))
    return(list(success = FALSE, subject_id = subject_id))
  }
}

# Process subjects (either sequentially or in parallel)
start_time <- Sys.time()

log_message("========================================")
log_message("FIRST PASS: Processing all subjects")
log_message("========================================")

if (opt$parallel && !opt$dry_run) {
  log_message(sprintf("Processing subjects in parallel with %d cores", opt$cores))
  
  # Create cluster
  cl <- makeCluster(min(opt$cores, total_subjects))
  
  # Export necessary variables and functions to the cluster
  clusterExport(cl, c("log_message", "SCRIPT_DIR", "opt", "fit_params", "data_params", "get_safe_data_dir"))
  
  # Process subjects in parallel
  results <- parLapply(cl, subject_indices, function(index) {
    fit_single_sub(index, subject_indices, subject_ids, opt, fit_params, data_params)
  })
  
  # Close cluster
  stopCluster(cl)
  
} else {
  # Process subjects sequentially
  results <- list()
  
  for (i in seq_along(subject_indices)) {
    index <- subject_indices[i]
    log_message(sprintf("Processing subject %d of %d", i, total_subjects))
    
    result <- fit_single_sub(index, subject_indices, subject_ids, opt, fit_params, data_params)
    results[[i]] <- result
  }
}

# Extract successful and failed indices
first_pass_success <- sapply(results, function(x) !is.null(x) && x$success)
first_pass_successful <- sum(first_pass_success)
first_pass_failed_indices <- subject_indices[!first_pass_success]
first_pass_failed <- length(first_pass_failed_indices)

log_message(sprintf("First pass complete: %d successful, %d failed", 
                    first_pass_successful, first_pass_failed))

# Retry failed subjects if retry is enabled
retry_successful <- 0
retry_failed <- 0
permanently_failed_indices <- c()

if (first_pass_failed > 0 && !opt$no_retry && !opt$dry_run) {
  log_message("")
  log_message("========================================")
  log_message(sprintf("RETRY PASS: Retrying %d failed subjects", first_pass_failed))
  log_message("========================================")
  
  # Optional: Add a brief delay to let system resources settle
  Sys.sleep(5)
  
  if (opt$parallel) {
    log_message(sprintf("Processing retry subjects in parallel with %d cores", opt$cores))
    
    # Create cluster
    cl <- makeCluster(min(opt$cores, first_pass_failed))
    
    # Export necessary variables and functions to the cluster
    clusterExport(cl, c("log_message", "SCRIPT_DIR", "opt", "fit_params", "data_params", "get_safe_data_dir"))
    
    # Process failed subjects in parallel
    retry_results <- parLapply(cl, first_pass_failed_indices, function(index) {
      fit_single_sub(index, subject_indices, subject_ids, opt, fit_params, data_params)
    })
    
    # Close cluster
    stopCluster(cl)
    
  } else {
    # Process failed subjects sequentially
    retry_results <- list()
    
    for (i in seq_along(first_pass_failed_indices)) {
      index <- first_pass_failed_indices[i]
      log_message(sprintf("Retrying subject %d of %d", i, first_pass_failed))
      
      result <- fit_single_sub(index, subject_indices, subject_ids, opt, fit_params, data_params)
      retry_results[[i]] <- result
    }
  }
  
  # Extract retry statistics
  retry_success <- sapply(retry_results, function(x) !is.null(x) && x$success)
  retry_successful <- sum(retry_success)
  retry_failed <- first_pass_failed - retry_successful
  permanently_failed_indices <- first_pass_failed_indices[!retry_success]
  
  log_message(sprintf("Retry pass complete: %d recovered, %d permanently failed", 
                      retry_successful, retry_failed))
  
  if (retry_failed > 0) {
    failed_subject_ids <- sapply(retry_results[!retry_success], function(x) x$subject_id)
    log_message(sprintf("Permanently failed subjects: %s", 
                        paste(failed_subject_ids, collapse = ", ")))
  }
}

# Checkpoint cleanup
if (!opt$keep_checkpoints && !opt$dry_run) {
  log_message("")
  log_message("========================================")
  log_message("CHECKPOINT CLEANUP")
  log_message("========================================")
  
  # Get all successful subject indices (first pass + retry)
  successful_indices <- c(
    subject_indices[first_pass_success],
    if (retry_successful > 0) first_pass_failed_indices[sapply(retry_results, function(x) !is.null(x) && x$success)] else c()
  )
  
  cleaned_count <- 0
  kept_count <- 0
  
  for (index in subject_indices) {
    subject_id <- subject_ids[index]
    if (is.na(subject_id) || subject_id == "") next
    
    checkpoint_path <- get_checkpoint_path(subject_id, opt)
    
    if (file.exists(checkpoint_path)) {
      if (index %in% successful_indices) {
        # Delete checkpoint for successful subject
        if (file.remove(checkpoint_path)) {
          cleaned_count <- cleaned_count + 1
          log_message(sprintf("Cleaned checkpoint for subject %s", subject_id))
        } else {
          log_message(sprintf("WARNING: Failed to delete checkpoint for subject %s", subject_id))
        }
      } else {
        # Keep checkpoint for failed subject
        kept_count <- kept_count + 1
        log_message(sprintf("Keeping checkpoint for failed subject %s at: %s", 
                            subject_id, checkpoint_path))
      }
    }
  }
  
  log_message(sprintf("Checkpoint cleanup complete: %d cleaned, %d kept for debugging", 
                      cleaned_count, kept_count))
}

# Combine batch results
total_successful <- first_pass_successful + retry_successful
if (!opt$dry_run && total_successful > 0) {
  log_message("")
  log_message("========================================")
  log_message("COMBINING BATCH DATA")
  log_message("========================================")
  
  # Build command arguments for combine script
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

log_message("")
log_message("========================================")
log_message("FINAL SUMMARY")
log_message("========================================")
log_message(sprintf("Total subjects processed: %d", total_subjects))
log_message(sprintf("First pass: %d successful, %d failed", first_pass_successful, first_pass_failed))
if (!opt$no_retry && first_pass_failed > 0 && !opt$dry_run) {
  log_message(sprintf("Retry pass: %d recovered, %d permanently failed", retry_successful, retry_failed))
}
log_message(sprintf("Final results: %d successful, %d failed", 
                    first_pass_successful + retry_successful, 
                    retry_failed))
log_message(sprintf("Total time: %.2f minutes", as.numeric(elapsed)))

cat("Script completed. Exiting explicitly.\n")