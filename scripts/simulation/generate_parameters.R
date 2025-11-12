#!/usr/bin/env Rscript

# Generate and save parameter sets
suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(optparse)
  library(here)
})

# Source helper functions
source(file.path(here::here(), "scripts", "parameter_recovery", "helper_functions_PR.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))

# Define command line options
option_list = list(
  make_option(c("-m", "--model"), type="character", help="Model name (e.g., ev, pvl)"),
  make_option(c("-t", "--task"), type="character", help="Task name (e.g., igt_mod, igt)"),
  make_option(c("-g", "--group"), type="character", default="sing", help="Group type"),
  make_option(c("--cohort"), type="character", default=NULL,
              help="Cohort identifier (e.g., ahrb, es)"),
  make_option(c("--session"), type="character", default=NULL,
              help="Session identifier [default: %default]"),
  make_option(c("-n", "--n_subjects"), type="integer", default=100),
  make_option(c("-d", "--method"), type="character", default="ibSPSepse"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="./Data/sim/params/"),
  make_option(c("-f", "--fit_file"), type="character", default=NULL, 
              help="Path to fit data (required for EPSE methods)"),
  make_option(c("-p", "--params"), type="character", default=NULL),
  make_option(c("--seed"), type="integer", default=NULL, help="Set seed. Default: random"),
  make_option(c("-e", "--exclude_file"), type="character", default=NULL, 
              help="Path to list of subs to filter out")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Options used:\n")
dput(opt)

# Set random seed for reproducibility
# If seed not provided, generate one
if (is.null(opt$seed)) {
  opt$seed <- sample.int(.Machine$integer.max, 1)
  message(sprintf("No seed provided. Using random seed: %d", opt$seed))
}
set.seed(opt$seed)

# Get directory structure
dirs <- setup_directories(opt$task)

# Set output directory
output_dir <- get_simulation_output_dir(opt$task, "parameters")
if (!is.null(opt$output_dir)) {
  output_dir <- opt$output_dir
}

# Make sure directory exists
ensure_dir_exists(output_dir)

# Source parameter generation functions
source(file.path(dirs$SIM_DIR, "param_gen.R"))

# Extract model name and create full model name
model_name <- opt$model
task <- opt$task
group_type <- opt$group
cohort <- opt$cohort
session <- opt$session
full_model_name <- paste(task, group_type, model_name, sep="_")

# Source required files to get model classes
source_required_files(dirs$SIM_DIR, opt$task)

# Initialize task and model objects
task_obj <- initialize_task(task, dirs$SIM_DIR)

if (grepl("batch", opt$group)){
  group_type = "sing"
} else {
  group_type = opt$group
}
model_obj <- initialize_model(model_name, task, task_obj, dirs$SIM_DIR, group_type)

# Handle model fit loading for empirical-based methods
if (opt$method %in% c("mbSPSepse", "sbSPSepse", "tSPSepse", "wpSPSepse", "hpsEPSE", "ibSPSepse")) {
  if (is.null(opt$fit_file)) {
    # Try to construct default path using BIDS-inspired filename
    fit_file <- file.path(dirs$DATA_DIR, "fits", "fit", opt$cohort, paste0("ses-", opt$session), 
                          generate_bids_filename(
                            prefix = NULL,
                            task = opt$task,
                            group = opt$group,
                            model = opt$model,
                            cohort = opt$cohort,
                            ses = opt$session,
                            additional_tags = list(
                              "type" = "fit",
                              "desc" = "output"
                            ),
                            ext = "rds"
                          ))
    
    if (!file.exists(fit_file)) {
      # Try without cohort/session as fallback
      fit_file <- file.path(dirs$DATA_DIR, "rds", "fit", 
                         paste0(task, "_", opt$group, "_", opt$model, "_fit_output.rds"))
      
      if (!file.exists(fit_file)) {
        stop(sprintf("Required model fit file not found: %s", fit_file))
      }
    }
  } else {
    fit_file <- opt$fit_file
    if (!file.exists(fit_file)) {
      stop(sprintf("Specified model fit file not found: %s", fit_file))
    }
  }
  
  model_fit <- tryCatch({
    fit <- readRDS(fit_file)
    # Detect batch files by checking if "batch" is in the fit_file path
    # This handles the case where GROUP_TYPE="sing" but the file is actually a batch file
    if (grepl("batch", fit_file) || grepl("batch", opt$group)) {
      fit$subjects <- names(fit)
      cat("Detected batch fit file. Added subjects field to fit object.\n")
    }
    fit
  }, error = function(e) {
    stop(sprintf("Error loading model fit file %s: %s", 
                 basename(fit_file), e$message))
  })
} else {
  model_fit <- NULL
}

# Filter data
subs_to_exclude <- NULL
if (!is.null(opt$exclude_file)){
  subs_to_exclude <- readLines(opt$exclude_file)
}

# Function to filter model_fit based on subs_to_exclude
filter_model_fit <- function(model_fit, subs_to_exclude) {
  if (is.null(subs_to_exclude) || length(subs_to_exclude) == 0 || is.null(model_fit)) {
    return(model_fit)
  }
  
  filtered_model_fit <- list()
  subs <- list()
  
  for (sublist_name in names(model_fit)) {
    if (sublist_name == "subjects") {
      next()
    }
    subid <- model_fit[[sublist_name]]$subid
    if (!subid %in% subs_to_exclude) {
      filtered_model_fit[[sublist_name]] <- model_fit[[sublist_name]]
      subs[[length(subs) + 1]] <- subid
    }
  }
  
  filtered_model_fit$subjects <- subs
  
  return(filtered_model_fit)
}

# Apply the filter function if model_fit is not NULL
if (!is.null(model_fit)) {
  model_fit <- filter_model_fit(model_fit, subs_to_exclude)
}

# Generate parameters using model object instead of config file
if (!is.null(model_fit) && "subjects" %in% names(model_fit)) {
  names(model_fit) <- c(unlist(model_fit$subjects), "subjects")
}

# Determine maximum number of subjects from model fit
if (!is.null(model_fit) && !is.null(model_fit$all_params)) {
  highest_num_subs <- max(as.numeric(gsub(".*\\[(\\d+)\\].*", "\\1", 
                                         model_fit$all_params[grepl("\\[\\d+\\]", model_fit$all_params)])), 
                         na.rm = TRUE)
  n_subs_to_generate <- min(opt$n_subjects, highest_num_subs)
} else {
  n_subs_to_generate <- opt$n_subjects
}

params <- generate_parameters(
  model = model_obj,
  method = opt$method,
  model_fit = model_fit,
  n_subjects = n_subs_to_generate
)

# Save parameters using BIDS-inspired filename
filename <- file.path(
  output_dir,
  generate_bids_filename(
    prefix = NULL,
    task = task,
    group = group_type,
    model = model_name,
    cohort = cohort,
    ses = session,
    additional_tags = list(
      "type" = "params",
      "desc" = opt$method,
      "n" = n_subs_to_generate
    ),
    ext = "rds"
  )
)

saveRDS(params, filename)
cat("Parameters saved to:", "\n")
cat(filename, "\n")
