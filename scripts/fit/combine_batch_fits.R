#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(here)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-k", "--task"), type="character", help="Task name (e.g., igt_mod)"),
  make_option(c("-m", "--model"), type="character", help="Model name (e.g., ev_ddm)"),
  make_option(c("-g", "--group"), type="character", default="sing", help="Group identifier (e.g., sing, batch_001)"),
  make_option(c("-s", "--source"), type="character", help="Data source (e.g., adb, es)"),
  make_option(c("--ses"), type="character", default=NULL, help="Session identifier (optional)"),
  make_option(c("-t", "--type"), type="character", default="fit", help="Model type (fit, postpc, prepc)"),
  make_option(c("-c", "--combine_batches"), type="logical", default=FALSE,
              action="store_true", help="Combine all existing batch files"),
  make_option(c("-b", "--combine_batches_only"), type="logical", default=FALSE,
              action="store_true", help="Only combine batch files (ignore individual files)"),
  make_option(c("-u", "--update_batch"), type="character", default="",
              help="Update specific batch (e.g., batch_001)"),
  make_option(c("-d", "--delete_source"), type="logical", default=FALSE,
              action="store_true", help="Delete source files after successful combination"),
  make_option(c("-n", "--dry_run"), type="logical", default=FALSE,
              action="store_true", help="Show what would be done without doing it"),
  make_option(c("--warn_checkpoints"), type="logical", default=TRUE,
              action="store_true", help="Warn about checkpoint files found (default: TRUE)")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Validate required arguments and options
if (is.null(opt$task) || is.null(opt$model) || is.null(opt$source)) {
  stop("Task, model name, and source are required")
}

# More explicit type checking and comparison
if (isTRUE(opt$combine_batches_only) && nchar(opt$update_batch) > 0) {
  stop("Cannot use --combine_batches_only with --update_batch")
}

# Load helper functions for directory structure
source(file.path(here::here(), "scripts", "helpers", "helper_dirs.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))

# Set up task-specific directories
rds_dir <- get_fits_output_dir(opt$task, opt$type, opt$source, opt$ses)
if (!dir.exists(rds_dir)) {
  dir.create(rds_dir, recursive = TRUE)
  cat("Created directory: ", rds_dir, "\n")
}

# Function to safely read RDS files
safe_read_rds <- function(file) {
  tryCatch({
    readRDS(file)
  }, error = function(e) {
    stop("Error reading file ", basename(file), ": ", e$message)
  })
}

# Function to get next batch number using BIDS naming
get_next_batch_num <- function(rds_dir, task, model, source=NULL, ses=NULL) {
  # Create pattern to match BIDS format batch files
  ses_part <- if (!is.null(ses)) sprintf("ses-%s_", ses) else ""
  pattern <- sprintf("task-%s_cohort-%s_%sgroup-batch_[0-9]{3}_model-%s_.*\\.rds$", task, source, ses_part, model)
  existing_files <- list.files(rds_dir, pattern = pattern)
  
  if (length(existing_files) == 0) {
    return(1)
  }
  
  # Extract batch numbers from BIDS filenames
  batch_nums <- as.numeric(str_extract(existing_files, "(?<=group-batch_)\\d{3}"))
  return(max(batch_nums) + 1)
}

# Function to find checkpoint files
find_checkpoint_files <- function(rds_dir, task, model, group="sing", type="fit", source=NULL, ses=NULL) {
  # Pattern matches checkpoint files
  ses_part <- if (!is.null(ses)) sprintf("ses-%s_", ses) else ""
  pattern <- sprintf("^task-%s_cohort-%s_%sgroup-%s_model-%s_sub-\\d+.*_checkpoint\\.rds$", 
                     task, source, ses_part, group, model)
  files <- list.files(rds_dir, pattern = pattern, full.names = TRUE)
  return(files)
}

# Function to find individual RDS files using BIDS naming (EXCLUDING checkpoints)
find_individual_files <- function(rds_dir, task, model, group="sing", type="fit", source=NULL, ses=NULL) {
  # Pattern matches BIDS format: task-{task}_cohort-{source}_[ses-{ses}_]group-{group}_model-{model}_sub-XXX_*
  # but NOT files ending with _checkpoint.rds
  ses_part <- if (!is.null(ses)) sprintf("ses-%s_", ses) else ""
  pattern <- sprintf("^task-%s_cohort-%s_%sgroup-%s_model-%s_sub-\\d+.*\\.rds$", 
                     task, source, ses_part, group, model)
  all_files <- list.files(rds_dir, pattern = pattern, full.names = TRUE)
  
  # Filter out checkpoint files
  checkpoint_pattern <- "_checkpoint\\.rds$"
  files <- all_files[!grepl(checkpoint_pattern, all_files)]
  
  return(files)
}

# Function to find batch files using BIDS naming
find_batch_files <- function(rds_dir, task, model, type="fit", source=NULL, ses=NULL) {
  # Pattern matches BIDS format: task-{task}_cohort-{source}_[ses-{ses}_]group-batch_XXX_model-{model}_*
  ses_part <- if (!is.null(ses)) sprintf("ses-%s_", ses) else ""
  pattern <- sprintf("task-%s_cohort-%s_%sgroup-batch_[0-9]{3}_model-%s.*\\.rds$", task, source, ses_part, model)
  files <- list.files(rds_dir, pattern = pattern, full.names = TRUE)
  return(files)
}

# Function to check for orphaned checkpoints (checkpoints with corresponding successful fits)
check_orphaned_checkpoints <- function(checkpoint_files, individual_files) {
  orphaned <- c()
  
  for (checkpoint_file in checkpoint_files) {
    # Remove _checkpoint.rds suffix to get the base filename
    base_name <- sub("_checkpoint\\.rds$", ".rds", basename(checkpoint_file))
    
    # Check if a successful fit exists
    matching_fit <- grep(base_name, basename(individual_files), fixed = TRUE, value = TRUE)
    
    if (length(matching_fit) > 0) {
      orphaned <- c(orphaned, checkpoint_file)
    }
  }
  
  return(orphaned)
}

# Function to load individual RDS files using BIDS naming
load_individual_files <- function(files) {
  if (length(files) == 0) {
    return(NULL)
  }
  
  result <- list()
  for (file in files) {
    # Extract subject number from BIDS format filename
    index <- str_extract(basename(file), "(?<=sub-)(\\d+)(?=_)")
    cat("Loading individual file for index", index, "\n")
    fit_data <- safe_read_rds(file)
    result[[index]] <- fit_data
  }
  
  return(result)
}

# Function to load batch files
load_batch_files <- function(files) {
  if (length(files) == 0) {
    return(NULL)
  }
  
  all_data <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    cat("Loading batch file:", basename(files[i]), "\n")
    all_data[[i]] <- safe_read_rds(files[i])
  }
  
  result <- do.call(c, all_data)
  return(result)
}

# Function to safely delete files
safe_delete_files <- function(files) {
  deleted_files <- c()
  failed_files <- c()
  
  for (file in files) {
    tryCatch({
      if (file.remove(file)) {
        deleted_files <- c(deleted_files, basename(file))
      } else {
        failed_files <- c(failed_files, basename(file))
      }
    }, error = function(e) {
      failed_files <<- c(failed_files, basename(file))
    })
  }
  
  list(deleted = deleted_files, failed = failed_files)
}

# Check for checkpoint files and warn if found
checkpoint_files <- find_checkpoint_files(rds_dir, opt$task, opt$model, opt$group, opt$type, opt$source, opt$ses)
if (length(checkpoint_files) > 0 && opt$warn_checkpoints) {
  cat("\n")
  cat("========================================\n")
  cat("WARNING: Checkpoint files detected\n")
  cat("========================================\n")
  cat(sprintf("Found %d checkpoint files in the directory.\n", length(checkpoint_files)))
  cat("These are temporary files from incomplete fits.\n")
  cat("They will NOT be included in the combined output.\n\n")
  
  # Check for orphaned checkpoints
  individual_files_all <- list.files(rds_dir, pattern = sprintf("^task-%s.*_sub-\\d+.*\\.rds$", opt$task), 
                                     full.names = TRUE)
  individual_files_all <- individual_files_all[!grepl("_checkpoint\\.rds$", individual_files_all)]
  orphaned <- check_orphaned_checkpoints(checkpoint_files, individual_files_all)
  
  if (length(orphaned) > 0) {
    cat(sprintf("Found %d orphaned checkpoints (successful fits exist):\n", length(orphaned)))
    for (file in orphaned) {
      cat(sprintf("  - %s\n", basename(file)))
    }
    cat("\nYou can safely delete these orphaned checkpoints.\n")
  }
  
  # Show checkpoints without successful fits (potential problem subjects)
  problem_checkpoints <- setdiff(checkpoint_files, orphaned)
  if (length(problem_checkpoints) > 0) {
    cat(sprintf("\nFound %d checkpoints without successful fits (may need retry):\n", length(problem_checkpoints)))
    for (file in problem_checkpoints) {
      cat(sprintf("  - %s\n", basename(file)))
    }
  }
  cat("========================================\n\n")
}

# Dry run output
if (opt$dry_run) {
  cat("\nDRY RUN - showing what would be processed:\n")
  cat("=========================================\n")
  cat("Parameters:\n")
  cat("  Task:", opt$task, "\n")
  cat("  Model:", opt$model, "\n")
  cat("  Source:", opt$source, "\n")
  if (!is.null(opt$ses)) cat("  Session:", opt$ses, "\n")
  cat("  Combine batches:", opt$combine_batches, "\n")
  cat("  Combine only:", opt$combine_batches_only, "\n")
  cat("  Update batch:", if (is.null(opt$update_batch) || opt$update_batch == "") "None" else opt$update_batch, "\n")
  cat("  Delete source:", opt$delete_source, "\n\n")
  
  if (!opt$combine_batches_only) {
    individual_files <- find_individual_files(rds_dir, opt$task, opt$model, opt$group, opt$type, opt$source, opt$ses)
    cat("Individual files found (", length(individual_files), "):\n")
    if (length(individual_files) > 0) {
      for (file in individual_files) {
        file_size <- file.size(file) / 1024 / 1024  # Convert to MB
        cat(sprintf("  - %s (%.2f MB)\n", basename(file), file_size))
      }
    } else {
      cat("  None\n")
    }
    cat("\n")
  }
  
  if (opt$combine_batches) {
    batch_files <- find_batch_files(rds_dir, opt$task, opt$model, opt$type, opt$source, opt$ses)
    cat("Batch files found (", length(batch_files), "):\n")
    if (length(batch_files) > 0) {
      for (file in batch_files) {
        file_size <- file.size(file) / 1024 / 1024  # Convert to MB
        cat(sprintf("  - %s (%.2f MB)\n", basename(file), file_size))
      }
    } else {
      cat("  None\n")
    }
    cat("\n")
  }
  
  # Show checkpoint files that will be excluded
  if (length(checkpoint_files) > 0) {
    cat("Checkpoint files found (WILL BE EXCLUDED) (", length(checkpoint_files), "):\n")
    for (file in checkpoint_files) {
      file_size <- file.size(file) / 1024 / 1024  # Convert to MB
      cat(sprintf("  - %s (%.2f MB)\n", basename(file), file_size))
    }
    cat("\n")
  }
  
  if (opt$update_batch != "") {
    # Use proper BIDS naming convention for the update file
    ses_part <- if (!is.null(opt$ses)) sprintf("ses-%s_", opt$ses) else ""
    update_filename <- sprintf("task-%s_cohort-%s_%sgroup-%s_model-%s_type-%s_desc-output.rds", 
                               opt$task, opt$source, ses_part, opt$update_batch, opt$model, opt$type)
    update_file <- file.path(rds_dir, update_filename)
    if (!file.exists(update_file)) {
      cat("WARNING: Update batch file does not exist:", basename(update_file), "\n\n")
    }
  }
  
  if (opt$update_batch != "") {
    # Use proper BIDS naming convention for the update file
    ses_part <- if (!is.null(opt$ses)) sprintf("ses-%s_", opt$ses) else ""
    output_file <- sprintf("task-%s_cohort-%s_%sgroup-%s_model-%s_type-%s_desc-output.rds", 
                           opt$task, opt$source, ses_part, opt$update_batch, opt$model, opt$type)
  } else {
    next_batch <- get_next_batch_num(rds_dir, opt$task, opt$model, opt$source, opt$ses)
    batch_name <- sprintf("batch_%03d", next_batch)
    ses_part <- if (!is.null(opt$ses)) sprintf("ses-%s_", opt$ses) else ""
    output_file <- sprintf("task-%s_cohort-%s_%sgroup-%s_model-%s_type-%s_desc-output.rds", 
                           opt$task, opt$source, ses_part, batch_name, opt$model, opt$type)
  }
  cat("Would save combined output to:\n")
  cat("  ", output_file, "\n")
  
  if (opt$delete_source) {
    # Collect files that would be deleted
    files_to_delete <- c()
    
    if (!opt$combine_batches_only) {
      individual_files <- find_individual_files(rds_dir, opt$task, opt$model, opt$group, opt$type, opt$source, opt$ses)
      files_to_delete <- c(files_to_delete, individual_files)
    }
    
    if (opt$combine_batches) {
      batch_files <- find_batch_files(rds_dir, opt$task, opt$model, opt$type, opt$source, opt$ses)
      files_to_delete <- c(files_to_delete, batch_files)
    }
    
    # Display files that would be deleted
    cat("Files that would be deleted:\n")
    if (length(files_to_delete) > 0) {
      for (file in files_to_delete) {
        file_size <- file.size(file) / 1024 / 1024  # Convert to MB
        cat(sprintf("  - %s (%.2f MB)\n", basename(file), file_size))
      }
      cat("Total files to be deleted:", length(files_to_delete), "\n")
    } else {
      cat("  None\n")
    }
    cat("\n")
  }
  
  quit(status = 0)
}

# Prepare to combine files
result_list <- list()
files_to_delete <- c()

# Load individual files (checkpoints are automatically excluded)
if (!opt$combine_batches_only) {
  individual_files <- find_individual_files(rds_dir, opt$task, opt$model, opt$group, opt$type, opt$source, opt$ses)
  individual_data <- load_individual_files(individual_files)
  if (!is.null(individual_data)) {
    result_list <- c(result_list, individual_data)
    cat("Loaded", length(individual_data), "individual files\n")
    
    # Track individual files for potential deletion
    if (opt$delete_source) {
      files_to_delete <- c(files_to_delete, individual_files)
    }
  } else {
    cat("No individual files found\n")
  }
}

# Load batch files
if (opt$combine_batches) {
  batch_files <- find_batch_files(rds_dir, opt$task, opt$model, opt$type, opt$source, opt$ses)
  batch_data <- load_batch_files(batch_files)
  if (!is.null(batch_data)) {
    result_list <- c(result_list, batch_data)
    cat("Added data from batch files\n")
    
    # Track batch files for potential deletion
    if (opt$delete_source) {
      files_to_delete <- c(files_to_delete, batch_files)
    }
  } else {
    cat("No batch files found\n")
  }
}

# Check if we have any data to process
if (length(result_list) == 0) {
  stop("No files found to process")
}

# Determine output file using BIDS naming
if (opt$update_batch != "") {
  # Create BIDS-compliant filename for the updated batch
  batch_name <- opt$update_batch
  output_filename <- generate_bids_filename(
    prefix = NULL,
    task = opt$task,
    cohort = opt$source,
    ses = opt$ses,
    group = batch_name,
    model = opt$model,
    additional_tags = list("type" = opt$type, "desc" = "output"),
    ext = "rds"
  )
  output_file <- file.path(rds_dir, output_filename)
  
  if (!file.exists(output_file)) {
    warning("Update batch file does not exist: ", basename(output_file))
  }
} else {
  # Create a new batch with incremented number
  batch_num <- get_next_batch_num(rds_dir, opt$task, opt$model, opt$source, opt$ses)
  batch_name <- sprintf("batch_%03d", batch_num)
  
  output_filename <- generate_bids_filename(
    prefix = NULL,
    task = opt$task,
    cohort = opt$source,
    ses = opt$ses,
    group = batch_name,
    model = opt$model,
    additional_tags = list("type" = opt$type, "desc" = "output"),
    ext = "rds"
  )
  output_file <- file.path(rds_dir, output_filename)
}

# Save combined data
saveRDS(result_list, output_file)
cat("Saved combined data to:", output_file, "\n")
cat("Total subjects included:", length(result_list), "\n")

# Delete source files if option is enabled
if (opt$delete_source && length(files_to_delete) > 0) {
  cat("\nAttempting to delete source files...\n")
  deletion_result <- safe_delete_files(files_to_delete)
  
  if (length(deletion_result$deleted) > 0) {
    cat("Successfully deleted files:\n")
    cat(paste("  -", deletion_result$deleted, collapse = "\n"), "\n")
  }
  
  if (length(deletion_result$failed) > 0) {
    warning("Failed to delete the following files:\n", 
            paste("  -", deletion_result$failed, collapse = "\n"))
  }
}