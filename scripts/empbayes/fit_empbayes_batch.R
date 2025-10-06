#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(parallel)
  library(dplyr)
  library(cmdstanr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="Model name"),
  make_option(c("-k", "--task"), type="character", default=NULL, 
              help="Task name"),
  make_option(c("-s", "--source"), type="character", default=NULL, 
              help="Data source (cohort)"),
  make_option(c("--ses"), type="character", default=NULL, 
              help="Session identifier (optional)"),
  make_option(c("--model_status"), type="character", default=NULL,
              help="Model status (canonical/experimental/working). Default: auto-detect"),
  make_option(c("-a", "--subjects"), type="character", default=NULL, 
              help="Subject indices (e.g., '40-45' or '40,42,45' or 'all')"),
  make_option(c("--n_subs"), type="integer", default=NULL, 
              help="Fit first N eligible subjects after filtering"),
  make_option(c("--n_trials"), type="integer", default=120, 
              help="Number of trials"),
  make_option(c("--RTbound_min_ms"), type="integer", default=50, 
              help="RT minimum bound in milliseconds"),
  make_option(c("--RTbound_max_ms"), type="integer", default=2500, 
              help="RT maximum bound in milliseconds"),
  make_option(c("--rt_method"), type="character", default="remove", 
              help="RT method"),
  make_option(c("--n_warmup"), type="integer", default=3000, 
              help="Number of warmup iterations"),
  make_option(c("--n_iter"), type="integer", default=15000, 
              help="Number of iterations"),
  make_option(c("--n_chains"), type="integer", default=4, 
              help="Number of chains"),
  make_option(c("--adapt_delta"), type="double", default=0.95, 
              help="Adapt delta"),
  make_option(c("--max_treedepth"), type="integer", default=12, 
              help="Max tree depth"),
  make_option(c("--check_iter"), type="integer", default=1000, 
              help="Iteration interval for checkpoint runs"),
  make_option(c("--seed"), type="integer", default=29518, 
              help="Random seed"),
  make_option(c("--priors_file"), type="character", default=NULL, 
              help="Path to priors CSV file (optional)"),
  make_option(c("-p", "--parallel"), action="store_true", default=FALSE, 
              help="Run subjects in parallel"),
  make_option(c("-c", "--cores"), type="integer", default=4, 
              help="Number of cores for parallel processing"),
  make_option(c("--dry_run"), action="store_true", default=FALSE, 
              help="Perform a dry run")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$model) || is.null(opt$task) || is.null(opt$source)) {
  stop("Model, task, and source are required. Use -m, -k, and -s options.")
}

if (is.null(opt$subjects) && is.null(opt$n_subs)) {
  stop("Either --subjects or --n_subs must be specified.")
}

# Load helper functions
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))

# Set random seed
set.seed(opt$seed)

# Create log directory
log_dir <- file.path(get_proj_dir(), "log_files", opt$task, "empbayes")
ensure_dir_exists(log_dir)

# Create log file
log_filename <- generate_bids_filename(
  prefix = "empbayes_batch",
  cohort = opt$source,
  ses = opt$ses,
  task = opt$task,
  model = opt$model,
  additional_tags = list("date" = format(Sys.time(), "%Y%m%d-%H%M%S")),
  ext = "log"
)
log_file <- file.path(log_dir, log_filename)

# Function to log messages
log_message <- function(message) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  log_entry <- sprintf("%s %s", timestamp, message)
  cat(log_entry, "\n")
  cat(log_entry, "\n", file = log_file, append = TRUE)
}

log_message("Starting Empirical Bayes batch fitting")
log_message(sprintf("Model: %s, Task: %s, Source: %s", opt$model, opt$task, opt$source))
if (!is.null(opt$ses)) log_message(sprintf("Session: %s", opt$ses))

# Function to expand array spec
expand_array_spec <- function(spec) {
  if (tolower(spec) == "all") {
    return("all")
  }
  
  result <- integer(0)
  ranges <- strsplit(spec, ",")[[1]]
  
  for (range in ranges) {
    if (grepl("-", range)) {
      range_parts <- as.integer(strsplit(range, "-")[[1]])
      result <- c(result, range_parts[1]:range_parts[2])
    } else {
      result <- c(result, as.integer(range))
    }
  }
  return(result)
}

# Determine priors file path
if (is.null(opt$priors_file)) {
  priors_dir <- file.path(get_proj_dir(), "Outputs", opt$task, "empbayes", "priors")
  priors_filename <- generate_bids_filename(
    prefix = NULL,
    task = opt$task,
    cohort = opt$source,
    ses = opt$ses,
    model = opt$model,
    additional_tags = list("priors" = ""),
    ext = "csv"
  )
  opt$priors_file <- file.path(priors_dir, priors_filename)
}

log_message(sprintf("Priors file: %s", opt$priors_file))

# Load priors
if (!file.exists(opt$priors_file) && !opt$dry_run) {
  stop("Priors file not found: ", opt$priors_file, "\nRun generate_empbayes_priors.R first.")
}

if (!opt$dry_run) {
  priors_df <- read.csv(opt$priors_file)
  log_message("Loaded priors:")
  print(priors_df)
  
  # Convert to list format expected by fit_and_save_model
  priors_list <- lapply(1:nrow(priors_df), function(i) {
    c(priors_df$pr_mu[i], priors_df$pr_sigma[i])
  })
  names(priors_list) <- priors_df$parameter
  
} else {
  log_message("Dry run: Would load priors from file")
  priors_list <- NULL
}

# Load data and determine subjects to fit
if (!opt$dry_run) {
  log_message("Loading data...")
  all_data <- load_data(opt$task, opt$source, opt$ses)
  all_subjects <- unique(all_data$subjID)
  
  # Apply same data quality filters as subject selection
  log_message("Applying data quality filters...")
  eligible_subjects <- c()
  
  for (subj in all_subjects) {
    subj_data <- all_data[all_data$subjID == subj, ]
    
    if (nrow(subj_data) < opt$n_trials) next
    
    if (opt$rt_method == "remove") {
      valid_trials <- subj_data$RT >= (opt$RTbound_min_ms / 1000) & 
                      subj_data$RT <= (opt$RTbound_max_ms / 1000)
      if (sum(valid_trials) < opt$n_trials) next
    }
    
    eligible_subjects <- c(eligible_subjects, subj)
  }
  
  log_message(sprintf("Eligible subjects: %d", length(eligible_subjects)))
  
  # Determine which subjects to fit
  if (!is.null(opt$subjects)) {
    subjects_spec <- expand_array_spec(opt$subjects)
    if (is.character(subjects_spec) && subjects_spec == "all") {
      subjects_to_fit <- eligible_subjects
    } else {
      subjects_to_fit <- eligible_subjects[subjects_spec]
      subjects_to_fit <- subjects_to_fit[!is.na(subjects_to_fit)]
    }
  } else {
    # Use first n_subs eligible subjects
    subjects_to_fit <- eligible_subjects[1:min(opt$n_subs, length(eligible_subjects))]
  }
  
  log_message(sprintf("Subjects to fit: %d", length(subjects_to_fit)))
  
} else {
  log_message("Dry run: Would load data and filter subjects")
  subjects_to_fit <- paste0("sub", 1:10)  # Dummy
}

# Function to fit single subject
fit_single_subject <- function(subid, opt, priors_list) {
  log_message(sprintf("Processing subject: %s", subid))
  
  # Load subject data
  all_data <- load_data(opt$task, opt$source, opt$ses)
  subject_data <- all_data[all_data$subjID == subid, ]
  
  if (nrow(subject_data) == 0) {
    log_message(sprintf("ERROR: No data for subject %s", subid))
    return(FALSE)
  }
  
  # Extract data
  model_defaults <- get_model_defaults(opt$task)
  full_model_name <- paste(opt$task, "sing", opt$model, sep="_")
  data_to_extract <- model_defaults[[full_model_name]]$data
  model_params <- model_defaults[[full_model_name]]$params
  
  data_list <- extract_sample_data(
    subject_data, 
    data_to_extract, 
    task = opt$task,
    n_trials = opt$n_trials, 
    RTbound_min_ms = opt$RTbound_min_ms,
    RTbound_max_ms = opt$RTbound_max_ms,
    RTbound_reject_min_ms = opt$RTbound_min_ms + 20,
    RTbound_reject_max_ms = opt$RTbound_max_ms, 
    rt_method = opt$rt_method,
    minrt_ep_ms = 0,
    SID = TRUE
  )
  
  # Handle entropy data
  if (sum(grepl("entropy", names(data_list))) > 0) {
    if ("shannon_entropy" %in% names(data_list)) {
      data_list$shannon_entropy <- as.vector(data_list$shannon_entropy)
    }
    if ("transition_entropy" %in% names(data_list)) {
      data_list$transition_entropy <- as.vector(data_list$transition_entropy)
    }
    if ("ngram_entropy" %in% names(data_list)) {
      data_list$ngram_entropy <- as.vector(data_list$ngram_entropy)
    }
    if ("conditional_entropy" %in% names(data_list)) {
      data_list$conditional_entropy <- as.vector(data_list$conditional_entropy)
    }
  }
  
  # Set output directory to empbayes/individual
  output_dir <- get_empbayes_output_dir(opt$task, "individual", opt$source)
  
  # Fit model with empirical Bayes priors
  tryCatch({
    fit <- fit_and_save_model(
      task = opt$task,
      cohort = opt$source,
      ses = opt$ses,
      group_type = "emp",
      model_name = opt$model,
      model_type = "fit",
      data_list = data_list,
      n_subs = 1,
      n_trials = opt$n_trials,
      n_warmup = opt$n_warmup,
      n_iter = opt$n_iter,
      n_chains = opt$n_chains,
      adapt_delta = opt$adapt_delta,
      max_treedepth = opt$max_treedepth,
      model_params = model_params,
      dry_run = FALSE,
      checkpoint_interval = opt$check_iter,
      output_dir = output_dir,
      emp_bayes = TRUE,
      informative_priors = priors_list,
      subid = subid,
      index = NULL,
      init_params = NULL,
      model_status = opt$model_status,
      cohort_sub_dir = FALSE
    )
    
    log_message(sprintf("Successfully fitted subject: %s", subid))
    return(TRUE)
    
  }, error = function(e) {
    log_message(sprintf("ERROR fitting subject %s: %s", subid, e$message))
    return(FALSE)
  })
}

# Fit subjects
start_time <- Sys.time()

if (opt$dry_run) {
  log_message("Dry run: Would fit subjects with empirical Bayes priors")
  log_message(sprintf("Would fit %d subjects", length(subjects_to_fit)))
  successful <- 0
  failed <- 0
} else {
  if (opt$parallel) {
    log_message(sprintf("Fitting subjects in parallel with %d cores", opt$cores))
    
    cl <- makeCluster(min(opt$cores, length(subjects_to_fit)))
    clusterExport(cl, c("log_message", "opt", "priors_list"))
    clusterEvalQ(cl, {
      library(here)
      source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
      source(file.path(here::here(), "scripts", "helpers", "helper_common.R"))
    })
    
    results <- parLapply(cl, subjects_to_fit, function(subid) {
      fit_single_subject(subid, opt, priors_list)
    })
    
    stopCluster(cl)
    
    successful <- sum(unlist(results))
    failed <- length(subjects_to_fit) - successful
    
  } else {
    log_message("Fitting subjects sequentially")
    successful <- 0
    failed <- 0
    
    for (i in seq_along(subjects_to_fit)) {
      subid <- subjects_to_fit[i]
      log_message(sprintf("Subject %d of %d", i, length(subjects_to_fit)))
      
      result <- fit_single_subject(subid, opt, priors_list)
      if (result) {
        successful <- successful + 1
      } else {
        failed <- failed + 1
      }
    }
  }
}

# Combine and move results
if (!opt$dry_run && successful > 0) {
  log_message("Combining individual fits...")
  
  # Find individual fit files
  individual_dir <- get_empbayes_output_dir(opt$task, "individual", opt$source)
  ses_part <- if (!is.null(opt$ses)) sprintf("ses-%s_", opt$ses) else ""
  pattern <- sprintf("^task-%s_cohort-%s_%sgroup-emp_model-%s_sub-.*\\.rds$", 
                     opt$task, opt$source, ses_part, opt$model)
  
  individual_files <- list.files(individual_dir, pattern = pattern, full.names = TRUE)
  log_message(sprintf("Found %d individual fit files", length(individual_files)))
  
  # Load and combine
  combined_fits <- list()
  for (file in individual_files) {
    fit_data <- readRDS(file)
    subid <- fit_data$subid
    combined_fits[[subid]] <- fit_data
  }
  
  # Save combined fits to final location
  final_output_dir <- get_fits_output_dir(opt$task, "fit", opt$source, opt$ses)
  ensure_dir_exists(final_output_dir)
  
  combined_filename <- generate_bids_filename(
    prefix = NULL,
    task = opt$task,
    cohort = opt$source,
    ses = opt$ses,
    group = "emp",
    model = opt$model,
    additional_tags = list("type" = "fit", "desc" = "output"),
    ext = "rds"
  )
  
  combined_file <- file.path(final_output_dir, combined_filename)
  saveRDS(combined_fits, combined_file)
  log_message(sprintf("Saved combined fits to: %s", combined_file))
  
  # Delete individual files
  log_message("Removing temporary individual fit files...")
  for (file in individual_files) {
    file.remove(file)
  }
  log_message("Cleanup complete")
}

# Log completion
end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")
log_message(sprintf("Batch fitting complete: %d successful, %d failed, %.2f minutes", 
                    successful, failed, as.numeric(elapsed)))
