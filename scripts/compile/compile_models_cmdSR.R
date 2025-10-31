#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(cmdstanr)
  library(digest)
})

#Sys.setenv("STAN_OPENCL" = "FALSE")
options(mc.cores = parallel::detectCores())

# Set up command line options
option_list <- list(
  make_option(c("-t", "--type"), type="character", default="all",
              help="Model type(s) to compile: fit, postpc, prepc, or all [default= %default]"),
  make_option(c("-k", "--task"), type="character", default=NULL,
              help="Task for model [default= %default]"),
  make_option(c("-o", "--model_status"), type="character", default=NULL,
              help="Required: validation status of model [default= %default]"),
  make_option(c("-m", "--models"), type="character", default="all",
              help="Primary string to match model names, or 'all' [default= %default]"),
  make_option(c("-s", "--secondary"), type="character", default=NULL,
              help="Comma-separated list of secondary strings to match model names (AND condition) [default= %default]"),
  make_option(c("-e", "--exclude"), type="character", default=NULL,
              help="Comma-separated list of strings to exclude from model names [default= %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="Print verbose output [default= %default]"),
  make_option(c("-n", "--dry-run"), action="store_true", default=FALSE,
              help="Show what would be compiled without actually compiling [default= %default]"),
  make_option(c("-y", "--yes"), action="store_true", default=FALSE,
              help="Automatically compile all matching models without prompting [default= %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Get task-specific directories
if (is.null(opt$model_status)) {
  stop("Validation status of model is required using the -o/--model_status option.")
}

# Load helper functions
source(file.path(here::here(), "scripts", "helpers", "helper_dirs.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))

# Set up directory paths
PROJ_DIR <- here::here()

if (is.null(opt$task)) {
  stop("Task parameter (-k) is required for compile_models_cmdSR.R")
} else {
  # Task-specific directories
  MODELS_DIR <- file.path(PROJ_DIR, "models", opt$task)
  MODELS_TXT_DIR <- file.path(MODELS_DIR, opt$model_status, "txt")
  MODELS_BIN_DIR <- file.path(MODELS_DIR, opt$model_status, "bin")
}

# Function to compile Stan model
compile_stan_model <- function(stan_file, bin_dir, output_filename, verbose = FALSE) {
  tryCatch({
    stan_filename <- basename(stan_file)
    
    # Set the exe_file to the output path
    exe_file <- file.path(bin_dir, output_filename)
    hash_file <- paste0(tools::file_path_sans_ext(exe_file), ".hash")
    
    # Calculate hash of the Stan file
    current_hash <- digest::digest(stan_file, file = TRUE)
    
    # Check if compilation is necessary
    need_compile <- TRUE
    if (file.exists(exe_file) && file.exists(hash_file)) {
      old_hash <- readLines(hash_file)
      if (current_hash == old_hash) {
        need_compile <- FALSE
        if (verbose) cat(paste0("No changes detected, skipping compilation: ", stan_filename, "\n"))
      }
    }
    
    if (need_compile) {
      if (verbose) cat(paste0("Compiling: ", stan_filename, "\n"))
      model <- cmdstan_model(stan_file, compile = TRUE,
                             force_recompile = TRUE, exe_file = exe_file,
                             quiet = !verbose)
      # Save the new hash
      writeLines(current_hash, hash_file)
      if (verbose) cat(paste0("Completed: ", output_filename, "\n\n"))
    }
  }, error = function(e) {
    cat(paste0("Error compiling ", stan_filename, ": ", e$message, "\n"))
    return(NULL)
  })
}

# Function to find matching models with support for comma-separated patterns
find_matching_models <- function(model_types, primary_pattern, secondary_pattern = NULL, exclude_pattern = NULL) {
  matching_models <- list()
  
  # Split patterns into lists if they contain commas
  secondary_patterns <- if (!is.null(secondary_pattern)) {
    strsplit(secondary_pattern, ",")[[1]]
  } else {
    NULL
  }
  
  exclude_patterns <- if (!is.null(exclude_pattern)) {
    strsplit(exclude_pattern, ",")[[1]]
  } else {
    NULL
  }
  
  for (type in model_types) {
    txt_dir <- file.path(MODELS_TXT_DIR, type)
    if (!dir.exists(txt_dir)) {
      warning(paste("Directory does not exist:", txt_dir))
      next
    }
    
    stan_files <- list.files(txt_dir, pattern = "\\.stan$", full.names = TRUE)
    
    if (primary_pattern != "all") {
      stan_files <- stan_files[grepl(primary_pattern, basename(stan_files), ignore.case = TRUE)]
    }
    
    # Apply all secondary patterns (AND condition)
    if (!is.null(secondary_patterns)) {
      for (pattern in secondary_patterns) {
        stan_files <- stan_files[grepl(pattern, basename(stan_files), ignore.case = TRUE)]
      }
    }
    
    # Apply all exclude patterns (OR condition)
    if (!is.null(exclude_patterns)) {
      exclude_mask <- sapply(basename(stan_files), function(x) {
        !any(sapply(exclude_patterns, function(pattern) grepl(pattern, x, ignore.case = TRUE)))
      })
      stan_files <- stan_files[exclude_mask]
    }
    
    matching_models[[type]] <- stan_files
  }
  
  return(matching_models)
}

# Main execution
main <- function() {
  # Determine model types to compile
  model_types <- if (opt$type == "all") c("fit", "postpc", "prepc") else strsplit(opt$type, ",")[[1]]
  
  # Find matching models
  matching_models <- find_matching_models(model_types, opt$models, opt$secondary, opt$exclude)
  
  # Check if any models were found
  total_models <- sum(sapply(matching_models, length))
  if (total_models == 0) {
    stop("No matching models found.")
  }
  
  # If multiple models found, show them and check if we should compile all
  if (total_models > 1 && opt$models != "all") {
    cat("Multiple matching models found:\n")
    for (type in names(matching_models)) {
      for (model in matching_models[[type]]) {
        cat(paste0("- ", type, ": ", basename(model), "\n"))
      }
    }
    if (!opt$yes && !opt[["dry-run"]]) {
      cat("To compile these models, re-run with the --yes option.\n")
      return(invisible(NULL))
    } else if (!opt$yes && opt[["dry-run"]]) {
      cat("To compile these models, would require running with the --yes option. Otherwise, it will exit unless you use a more precise model str.\n")
      return(invisible(NULL))
    }
  }
  
  # Compile models or show what would be compiled
  for (type in names(matching_models)) {
    for (stan_file in matching_models[[type]]) {
      type_bin_dir <- file.path(MODELS_BIN_DIR, type)
      
      # Extract base name without extension
      base_name <- tools::file_path_sans_ext(basename(stan_file))
      
      # Determine output filename
      output_filename <- base_name
      
      if (grepl("task-", base_name)) {
        # Already in BIDS format, check if it has a type tag
        if (!grepl("type-", base_name)) {
          # Add the type tag if not present
          components <- parse_bids_filename(base_name)
          output_filename <- generate_bids_filename(
            prefix = NULL,
            task = components$task,
            group = components$group,
            model = components$model,
            additional_tags = c(components[!names(components) %in% c("task", "group", "model")], 
                                list(type = type)),
            ext = NULL  # No extension yet
          )
        }
      } else {
        # Not in BIDS format, convert to BIDS
        parts <- strsplit(base_name, "_")[[1]]
        
        if (length(parts) >= 3) {
          task_name <- opt$task
          group_type <- parts[2]
          model_name <- parts[3]
        } else {
          # Default values if we can't parse
          task_name <- opt$task
          group_type <- "sing"
          model_name <- base_name
        }
        
        output_filename <- generate_bids_filename(
          prefix = NULL,
          task = task_name,
          group = group_type,
          model = model_name,
          additional_tags = list(type = type),
          ext = NULL  # No extension yet
        )
      }
      
      # Add the correct extension
      output_file <- paste0(output_filename, ".stan")
      
      # Ensure directory exists
      dir.create(type_bin_dir, showWarnings = FALSE, recursive = TRUE)
      
      if (isTRUE(opt[["dry-run"]])) {
        cat(paste("Would compile:", basename(stan_file), "to", file.path(type_bin_dir, output_file), "\n"))
      } else {
        compile_stan_model(stan_file, type_bin_dir, output_file, verbose = opt$verbose)
      }
    }
  }
}

# Run the main function
main()