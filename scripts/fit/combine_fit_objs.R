#!/usr/bin/env Rscript

# ---- Load required libraries ----
suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

# ----Function definitions----
initialize_combined_diagnostics <- function(max_treedepth) {
  list(
    subject_summaries = list(),
    flags = list(
      subjects_with_divergences = integer(0),
      subjects_max_treedepth = integer(0),
      subjects_low_accept_stat = integer(0)
    ),
    total_subjects = 0,
    max_treedepth = max_treedepth
  )
}

update_combined_diagnostics <- function(combined_diagnostics, new_diagnostics, subject_index) {
  combined_diagnostics$total_subjects <- combined_diagnostics$total_subjects + 1
  
  subject_summary <- list()
  
  diag_names <- c("treedepth__", "energy__", "accept_stat__", "stepsize__", "n_leapfrog__")
  
  for (diag in diag_names) {
    new_vals <- as.vector(new_diagnostics[,,diag])
    new_vals <- new_vals[!is.na(new_vals)]  # Remove NA values
    
    if (length(new_vals) > 0) {
      subject_summary[[diag]] <- list(
        mean = mean(new_vals),
        sd = sd(new_vals),
        min = min(new_vals),
        max = max(new_vals)
      )
    }
  }
  
  # Handle divergences separately
  divergences <- as.vector(new_diagnostics[,,"divergent__"])
  divergences <- divergences[!is.na(divergences)]
  if (sum(divergences) > 0) {
    subject_summary$divergent__ <- sum(divergences)
    combined_diagnostics$flags$subjects_with_divergences <- c(
      combined_diagnostics$flags$subjects_with_divergences, 
      subject_index
    )
  }
  
  # Store subject summary
  combined_diagnostics$subject_summaries[[as.character(subject_index)]] <- subject_summary
  
  # Update flags
  if (!is.null(subject_summary$treedepth__) && subject_summary$treedepth__$max >= combined_diagnostics$max_treedepth) {
    combined_diagnostics$flags$subjects_max_treedepth <- c(
      combined_diagnostics$flags$subjects_max_treedepth, 
      subject_index
    )
  }
  
  if (!is.null(subject_summary$accept_stat__)) {
    if (subject_summary$accept_stat__$mean < 0.2) {
      combined_diagnostics$flags$subjects_low_accept_stat <- c(
        combined_diagnostics$flags$subjects_low_accept_stat, 
        subject_index
      )
    }
  }
  
  combined_diagnostics
}

get_final_summary <- function(combined_diagnostics) {
  summary <- list()
  
  diag_names <- c("treedepth__", "energy__", "accept_stat__", "stepsize__", "n_leapfrog__")
  
  for (diag in diag_names) {
    all_vals <- sapply(combined_diagnostics$subject_summaries, function(x) x[[diag]]$mean)
    summary[[diag]] <- list(
      mean = mean(all_vals),
      sd = sd(all_vals),
      min = min(sapply(combined_diagnostics$subject_summaries, function(x) x[[diag]]$min)),
      max = max(sapply(combined_diagnostics$subject_summaries, function(x) x[[diag]]$max))
    )
  }
  
  summary$overall <- list(
    total_subjects = combined_diagnostics$total_subjects,
    subjects_with_divergences = length(combined_diagnostics$flags$subjects_with_divergences),
    total_divergences = sum(sapply(combined_diagnostics$subject_summaries, function(x) ifelse(is.null(x$divergent__), 0, x$divergent__)))
  )
  
  summary$flags <- combined_diagnostics$flags
  
  summary
}

parse_indices <- function(indices_str) {
  if (indices_str == "all") return(NULL)
  indices <- unlist(strsplit(indices_str, ","))
  result <- integer()
  for (index in indices) {
    if (grepl("-", index)) {
      range <- as.integer(unlist(strsplit(index, "-")))
      result <- c(result, seq(range[1], range[2]))
    } else {
      result <- c(result, as.integer(index))
    }
  }
  sort(unique(result))
}

extract_subject_number <- function(filename) {
  subject_pattern <- "sub-(\\d+)_"
  match <- regexec(subject_pattern, basename(filename))
  if (match[[1]][1] != -1) {
    return(as.integer(regmatches(basename(filename), match)[[1]][2]))
  } else {
    return(NA)
  }
}

combine_fits <- function(files) {
  combined_fit <- list()
  skipped_files <- character()
  
  for (i in seq_along(files)) {
    if (is.na(files[i])[[1]]) {
      cat("Subject file", i, "not found.\n")
      next
      }
    if (opt$verbose) cat("Processing file", i, "of", length(files), "\n")
    
    fit <- readRDS(files[i])
    subject_number <- extract_subject_number(files[i])
    
    # Validate consistency and init diagnostics
    if (i == 1) {
      combined_diagnostics <- initialize_combined_diagnostics(fit$max_treedepth)
      
      combined_fit$n_warmup <- fit$n_warmup
      combined_fit$n_iter <- fit$n_iter
      combined_fit$n_chains <- fit$n_chains
      combined_fit$adapt_delta <- fit$adapt_delta
      combined_fit$max_treedepth <- fit$max_treedepth
      combined_fit$tss <- fit$tss
      combined_fit$cmdstan_version <- fit$cmdstan_version
      combined_fit$model_name <- fit$model_name
      combined_fit$model_params <- fit$model_params
    } else {
      if (!identical(combined_fit$n_warmup, fit$n_warmup) ||
          !identical(combined_fit$n_iter, fit$n_iter) ||
          !identical(combined_fit$n_chains, fit$n_chains) ||
          !identical(combined_fit$adapt_delta, fit$adapt_delta) ||
          !identical(combined_fit$max_treedepth, fit$max_treedepth) ||
          !identical(combined_fit$tss, fit$tss) ||
          !identical(combined_fit$cmdstan_version, fit$cmdstan_version) ||
          !identical(combined_fit$model_name, fit$model_name) ||
          !identical(combined_fit$model_params, fit$model_params)) {
        warning("Inconsistent fit object found in file: ", files[i], ". Skipping this file.")
        skipped_files <- c(skipped_files, files[i])
        next  # Skip to the next iteration of the loop
      }
    }
    
    # Update parameter names with subject index
    new_var_names <- paste0(dimnames(fit$draws)[[3]], "[", subject_number, "]")
    dimnames(fit$draws)[[3]] <- new_var_names
    
    # Combine draws
    if (i == 1) {
      combined_fit$draws <- fit$draws
    } else {
      combined_fit$draws <- posterior::bind_draws(combined_fit$draws, fit$draws, along = "variable")
    }
    
    # Combine sampler diagnostics
    combined_diagnostics <- update_combined_diagnostics(combined_diagnostics, fit$sampler_diagnostics, subject_number)
    
    # Combine and update other elements
    if (i == 1) {
      combined_fit$all_params <- fit$all_params
      combined_fit$list_params <- fit$list_params
      combined_fit$n_runs <- fit$n_runs
      combined_fit$subid <- fit$subid
      combined_fit$index <- fit$index
    } else {
      combined_fit$all_params <- c(combined_fit$all_params, paste0(fit$all_params, "[", subject_number, "]"))
      combined_fit$list_params <- unique(c(combined_fit$list_params, fit$list_params))
      combined_fit$n_runs <- combined_fit$n_runs + fit$n_runs
      combined_fit$subid <- c(combined_fit$subid, fit$subid)
      combined_fit$index <- c(combined_fit$index, fit$index)
    }
  }
  
  # Update n_params
  combined_fit$n_params <- dim(combined_fit$draws)[3]
  
  # Recalculate summary stats
  combined_fit$summary_stats <- posterior::summarize_draws(combined_fit$draws)
  
  # Generate new params sample
  combined_fit$params <- unname(extract_params(combined_fit$all_params, n_subs = length(files), main_params_vec = combined_fit$model_params))
  
  # After processing all files, get the final summary
  combined_fit$sampler_diagnostics <- get_final_summary(combined_diagnostics)
  
  if (length(skipped_files) > 0) {
    cat("The following files were skipped due to inconsistencies:\n")
    cat(paste(skipped_files, collapse = "\n"), "\n")
    cat("Total files skipped:", length(skipped_files), "\n")
  }
  
  combined_fit$skipped_files <- skipped_files
  
  return(combined_fit)
}

# ----Parse command line arguments----
option_list = list(
  make_option(c("-k", "--task"), type="character", default=NULL, help="Task name", metavar="character"),
  make_option(c("-m", "--model"), type="character", default=NULL, help="Model name", metavar="character"),
  make_option(c("-t", "--model-type"), type="character", default="fit", help="Model type (fit, postpc, prepc)", metavar="character"),
  make_option(c("-s", "--source"), type="character", default=NULL, help="Data source (e.g., adb, es)", metavar="character"),
  make_option(c("--ses"), type="character", default=NULL, help="Session identifier (optional)", metavar="character"),
  make_option(c("-i", "--input-dir"), type="character", default=NULL, help="Input directory", metavar="character"),
  make_option(c("-o", "--output-dir"), type="character", default=NULL, help="Output directory", metavar="character"),
  make_option(c("-d", "--delete-individual"), action="store_true", default=FALSE, help="Delete individual subject outputs after combining"),
  make_option(c("-n", "--indices"), type="character", default="all", help="Indices of subjects to include (e.g., '1-100,102,105-200')", metavar="character"),
  make_option(c("-b", "--batch-size"), type="integer", default=0, help="Number of subjects to process in each batch", metavar="number"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print detailed progress information")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$task) || is.null(opt$model) || is.null(opt$source)) {
  stop("Task, model, and source must be specified")
}

# ----Set up directories----
# Load helper functions
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))

if (is.null(opt$input_dir)) {
  opt$input_dir <- get_fits_output_dir(opt$task, opt$type, opt$source, opt$ses)
}
if (is.null(opt$output_dir)) {
  opt$output_dir <- opt$input_dir
}

# Get list of files to process
ses_part <- if (!is.null(opt$ses)) sprintf("ses-%s_", opt$ses) else ""
file_pattern <- sprintf("^task-%s_cohort-%s_%sgroup-sing_model-%s_sub-\\d+.*\\.rds$", 
                       opt$task, opt$source, ses_part, opt$model)
all_files <- list.files(opt$input_dir, pattern = file_pattern, full.names = TRUE)

if (length(all_files) == 0) {
  stop("No matching files found in the input directory: ", opt$input_dir)
}

# Create a named vector of files, with names being the subject numbers
file_map <- setNames(all_files, sapply(all_files, extract_subject_number))

# Remove any NA entries (files that didn't match the pattern)
file_map <- file_map[!is.na(names(file_map))]

# Parse indices
indices <- parse_indices(opt$indices)

# Select files based on indices
if (!is.null(indices)) {
  files_to_process <- file_map[as.character(indices)]
  missing_indices <- indices[!indices %in% as.numeric(names(files_to_process))]
} else {
  files_to_process <- file_map
  missing_indices <- integer(0)
}

if (length(files_to_process) == 0) {
  stop("No files to process after applying indices")
}

# Report on missing files
if (length(missing_indices) > 0) {
  cat("Warning: The following indices do not have corresponding files:\n")
  cat(paste(missing_indices, collapse = ", "), "\n")
}

# ----Process files----
if (opt$batch_size > 0) {
  batches <- split(files_to_process, ceiling(seq_along(files_to_process)/opt$batch_size))
  combined_fits <- list()
  all_skipped_files <- character()
  
  for (i in seq_along(batches)) {
    if (opt$verbose) cat("Processing batch", i, "of", length(batches), "\n")
    batch_fit <- combine_fits(batches[[i]])
    combined_fits[[i]] <- batch_fit
    all_skipped_files <- c(all_skipped_files, batch_fit$skipped_files)
    
    # Save intermediate result with BIDS-style naming
    ses_str <- if (!is.null(opt$ses)) sprintf("ses-%s_", opt$ses) else ""
    batch_output_file <- file.path(opt$output_dir, 
                                 sprintf("task-%s_cohort-%s_%sgroup-batch_%03d_model-%s_%s_parallel_output.rds", 
                                        opt$task, opt$source, ses_str, i, opt$model, opt$model_type))
    saveRDS(batch_fit, batch_output_file)
  }
  
  # Combine batches
  if (opt$verbose) cat("Combining batches\n")
  final_fit <- combine_fits(lapply(combined_fits, function(fit) fit$draws))
  final_fit$skipped_files <- unique(all_skipped_files)
} else {
  if (opt$verbose) cat("Processing all files\n")
  final_fit <- combine_fits(files_to_process)
}

final_fit$missing_indices = missing_indices

# ----Save final combined fit----
# Generate BIDS-style filename for the output
ses_str <- if (!is.null(opt$ses)) sprintf("ses-%s_", opt$ses) else ""
output_filename <- sprintf("task-%s_cohort-%s_%sgroup-parallel_model-%s_type-%s_desc-output.rds", 
                          opt$task, opt$source, ses_str, opt$model, opt$model_type)
output_file <- file.path(opt$output_dir, output_filename)
saveRDS(final_fit, output_file)

if (opt$verbose) cat("Saved combined output to:", output_file, "\n")

# If verbose output is enabled, add:
if (opt$verbose) {
  cat("Total files processed:", length(files_to_process) - length(final_fit$skipped_files) - length(missing_indices), "\n")
  cat("Total files skipped:", length(final_fit$skipped_files), "\n")
  cat("Indices without corresponding files:", length(missing_indices), "\n")
}

# Delete individual files if specified
if (opt$delete_individual) {
  if (opt$verbose) cat("Deleting individual subject files\n")
  files_to_remove = files_to_process[!files_to_process %in% final_fit$skipped_files]
  file.remove(files_to_remove)
}

if (opt$verbose) cat("Processing complete\n")
