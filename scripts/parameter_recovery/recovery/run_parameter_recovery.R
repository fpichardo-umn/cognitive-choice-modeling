#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
})

# Source helper functions - this loads all the other helper functions
source(file.path(here::here(), "scripts", "parameter_recovery", "helper_functions_PR.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))

# Define command line options
option_list = list(
  make_option(c("-D", "--sim_data"), type="character", help="Simulated data file - simulated trial data and true parameter values"),
  make_option(c("-m", "--model"), type="character", help="Model name (e.g., ev, pvl)"),
  make_option(c("-k", "--task"), type="character", help="Task name (e.g., igt_mod)"),
  make_option(c("-g", "--group"), type="character", default="sing", help="Group type"),
  make_option(c("--cohort"), type="character", default=NULL,
              help="Cohort identifier (e.g., ahrb, es)"),
  make_option(c("--session"), type="character", default=NULL,
              help="Session identifier [default: %default]"),
  make_option(c("-n", "--indiv"), action="store_true", default=FALSE),
  make_option(c("-f", "--output_fit_dir"), type="character", default=NULL),
  make_option(c("-r", "--output_rec_dir"), type="character", default=NULL),
  make_option(c("-w", "--n_warmup"), type="integer", default=1000),
  make_option(c("-i", "--n_iter"), type="integer", default=2000),
  make_option(c("-c", "--n_chains"), type="integer", default=4),
  make_option(c("-s", "--seed"), type="integer", default=12345),
  make_option(c("-a", "--adapt_delta"), type="double", default=0.95),
  make_option(c("-d", "--max_treedepth"), type="integer", default=12),
  make_option(c("--check_iter"), type="integer", default=10000, help="Iteration interval for checkpoint runs. Default: 10000"),
  make_option(c("-R", "--render"), action="store_true", default=FALSE, help="Render RMD to HTML"),
  make_option(c("--n_trials"), type="integer", default=120, help="Number of trials"),
  make_option(c("--rt_method"), type="character", default="mark", help="RT method"),
  make_option(c("--RTbound_min_ms"), type="integer", default=50, help="RT min bound in milliseconds"),
  make_option(c("--RTbound_max_ms"), type="integer", default=4000, help="RT max bound in milliseconds"),
  make_option(c("--min_valid_rt_pct"), type="double", default=0.7, help="Minimum percent valid RT")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

dput(opt)

# Get directory structure
dirs <- setup_directories(opt$task)

# Source required files
source(file.path(dirs$PR_DIR, "recovery", "recovery.R"))

# Set output directories
output_fit_dir <- get_validation_output_dir(opt$task, "parameter_recovery", "fits")
if (!is.null(opt$output_fit_dir)) {
  output_fit_dir <- opt$output_fit_dir
}

output_rec_dir <- get_validation_output_dir(opt$task, "parameter_recovery", "analysis")
if (!is.null(opt$output_rec_dir)) {
  output_rec_dir <- opt$output_rec_dir
}

# Make sure directories exist
ensure_dir_exists(output_fit_dir)
ensure_dir_exists(output_rec_dir)

# Load simulated data - use new path if not specified
if (is.null(opt$sim_data)) {
  sim_data_file <- file.path(
    get_simulation_output_dir(opt$task, "data"),
    "rds",
    generate_bids_filename(
      prefix = NULL,
      task = opt$task,
      group = opt$group,
      model = opt$model,
      cohort = opt$cohort,
      ses = opt$session,
      additional_tags = list(
        "type" = "sim",
        "desc" = "data"
      ),
      ext = "rds"
    )
  )
  
  if (!file.exists(sim_data_file)) {
    stop("Simulated data file not found: ", sim_data_file)
  }
  
  sim_data <- readRDS(sim_data_file)
} else {
  sim_data <- readRDS(opt$sim_data)
}

# Extract model name and create full model name
model_name <- opt$model
task <- opt$task
group_type <- opt$group
cohort <- opt$cohort
session <- opt$session
full_model_name <- paste(task, group_type, model_name, sep="_")

# Get model defaults
model_defaults <- get_model_defaults(task)

# Get data and parameters for this model
# For batch groups, strip the batch number for lookup
model_key <- sub("batch_[0-9]*", "sing", full_model_name)
data_to_extract <- model_defaults[[model_key]]$data
model_params <- model_defaults[[model_key]]$params

# Extract data from simulation
stan_data <- extract_simulation_data(
  data = sim_data, 
  data_params = data_to_extract, 
  task = opt$task, 
  individual = opt$indiv,
  n_trials = opt$n_trials,
  
  # --- Pass all RT parameters ---
  RTbound_min_ms = opt$RTbound_min_ms, 
  RTbound_max_ms = opt$RTbound_max_ms,
  
  # Use the same "reject" logic as the real data scripts
  RTbound_reject_min_ms = opt$RTbound_min_ms + 20, 
  RTbound_reject_max_ms = opt$RTbound_max_ms, 
  
  rt_method = opt$rt_method, 
  min_valid_rt_pct = opt$min_valid_rt_pct
)

if ("subjects" %in% names(stan_data)){
  stan_data$subjects = NULL
}

# Determine fitting group
# For batch groups, we fit individual models but can analyze them as a group
is_batch <- grepl("batch", opt$group)
if (is_batch) {
  fit_group <- "sing"
  batch_id <- sub("batch_", "", opt$group)
  message("Using batch group: ", batch_id)
} else {
  fit_group <- opt$group
  batch_id <- NULL
}

# Define the combined fit file name using BIDS-inspired format
combined_fit_file <- file.path(
  output_fit_dir,
  generate_bids_filename(
    prefix = NULL,
    task = task,
    group = group_type,
    model = model_name,
    cohort = cohort,
    ses = session,
    additional_tags = list(
      "type" = "rec",
      "desc" = "fits"
    ),
    ext = "rds"
  )
)

# Initialize for storing results
fits <- list()
individual_fit_files <- c()

if (opt$indiv || is_batch) {
  # Individual model fits
  for (s_idx in 1:length(names(stan_data))) {
    sub <- names(stan_data)[s_idx]
    
    # For batch groups, include batch identifier in filename
    if (is_batch) {
      subid <- paste0(sub, "_batch", batch_id, "_sim")
    } else {
      subid <- paste0(sub, "_sim")
    }
    
    # Define individual fit file path with BIDS-inspired format
    subject_fit_file <- file.path(
      output_fit_dir,
      generate_bids_filename(
        prefix = NULL,
        task = task,
        group = fit_group,
        model = model_name,
        cohort = cohort,
        ses = session,
        additional_tags = list(
          "sub" = sub,
          "type" = "fit",
          "desc" = "output"
        ),
        ext = "rds"
      )
    )
    
    individual_fit_files <- c(individual_fit_files, subject_fit_file)
    
    message(sprintf("Fitting model for subject %s (%d of %d)", 
                   sub, s_idx, length(names(stan_data))))
    
    sub_stan_data <- stan_data[[sub]]
    sub_stan_data$task = NULL
    sub_stan_data$sid <- as.numeric(sub)
    
    # Fit model for this subject
    fit <- fit_and_save_model(
      task = opt$task,
      cohort = cohort,
      ses = session,
      group_type = fit_group,
      model_name = opt$model,
      model_type = "fit",
      data_list = sub_stan_data,
      n_subs = 1,
      n_trials = sub_stan_data$T,
      n_warmup = opt$n_warmup,
      n_iter = opt$n_iter,
      n_chains = opt$n_chains,
      adapt_delta = opt$adapt_delta,
      max_treedepth = opt$max_treedepth,
      model_params = model_params,
      output_dir = output_fit_dir, 
      index = sub, 
      checkpoint_interval = opt$check_iter,
      is_simulation = TRUE,
      cohort_sub_dir = FALSE
    )
    
    # Save individual fit
    saveRDS(fit, subject_fit_file)
    
    # Add to list of fits
    fits[[length(fits) + 1]] <- fit
  }
  
  # Save combined fits
  message("Saving combined fit file: ", combined_fit_file)
  saveRDS(fits, combined_fit_file)
  
  # Verify combined fit saved successfully
  if (file.exists(combined_fit_file)) {
    # Remove individual fit files
    message("Removing individual fit files")
    file.remove(individual_fit_files)
  } else {
    warning("Failed to save combined fit file. Keeping individual fit files.")
  }
} else {
  # Hierarchical model fit
  message("Fitting hierarchical model")
  
  # Define the hierarchical fit file with BIDS-inspired format
  hier_fit_file <- file.path(
    output_fit_dir,
    generate_bids_filename(
      prefix = NULL,
      task = task,
      group = "hier",
      model = model_name,
      cohort = cohort,
      ses = session,
      additional_tags = list(
        "type" = "rec",
        "desc" = "fit"
      ),
      ext = "rds"
    )
  )
  
  individual_fit_files <- c(individual_fit_files, hier_fit_file)
  
  # Fit hierarchical model
  fit <- fit_and_save_model(
    task = opt$task,
    cohort = cohort,
    ses = session,
    group_type = fit_group,
    model_name = opt$model,
    model_type = "fit",
    data_list = stan_data,
    n_subs = stan_data$N,
    n_trials = stan_data$T,
    n_warmup = opt$n_warmup,
    n_iter = opt$n_iter,
    n_chains = opt$n_chains,
    adapt_delta = opt$adapt_delta,
    max_treedepth = opt$max_treedepth,
    model_params = model_params,
    output_dir = output_fit_dir, 
    checkpoint_interval = opt$check_iter, 
    cohort_sub_dir = FALSE
  )
  
  # Save hierarchical fit
  saveRDS(fit, hier_fit_file)
  
  fits <- fit  # For consistency with individual case
}

# Extract true parameters for comparison
sub_params <- sim_data %>%
  as.data.frame() %>% 
  group_by(idx) %>% 
  select(all_of(model_params)) %>% 
  unique()

# Create CSV output directory if it doesn't exist
if (!dir.exists(output_rec_dir)) {
  dir.create(output_rec_dir, recursive = TRUE)
}

# Define recovery data output CSV with BIDS-inspired format
recovery_csv_file <- file.path(
  output_rec_dir,
  generate_bids_filename(
    prefix = NULL,
    task = task,
    group = group_type,
    model = model_name,
    cohort = cohort,
    ses = session,
    additional_tags = list(
      "type" = "rec",
      "desc" = "data"
    ),
    ext = "csv"
  )
)

# Extract and save recovery data
if (opt$indiv || is_batch) {
  # Individual fits or batch group
  sub_true_params <- list()
  sub_rec_params <- list()
  
  for (i in 1:length(fits)) {
    # Get corresponding subject from parameters
    # Note: this assumes order of subjects in fits matches order in sub_params
    sub_true_params[[length(sub_true_params) + 1]] <- sub_params[i,] %>% 
      ungroup() %>% 
      as.data.frame() %>% 
      select(-c("idx"))
    
    current_params = extract_recovery_parameters(
      fit = fits[[i]],
      model_type = "individual_fit",
      parameters = model_params,
      summary_fn = median
    )
    sub_rec_params[[length(sub_rec_params) + 1]] <- setNames(
      as.numeric(current_params$value),
      current_params$parameter
    )
  }
  
  # Convert to standard format
  true_params_df <- do.call(rbind, sub_true_params)
  recovered_params_df <- as.data.frame(do.call(rbind, sub_rec_params))
  
  # Create CSV with recovery data with improved error handling
  tryCatch({
    # Validate that data exists
    if (ncol(true_params_df) == 0 || nrow(true_params_df) == 0) {
      stop("Empty true parameters data frame")
    }
    if (ncol(recovered_params_df) == 0 || nrow(recovered_params_df) == 0) {
      stop("Empty recovered parameters data frame")
    }
    
    # Ensure output directory exists
    output_dir_path <- dirname(recovery_csv_file)
    if (!dir.exists(output_dir_path)) {
      dir.create(output_dir_path, recursive = TRUE)
    }
    
    # Analyze recovery using the new approach and save directly to CSV
    recovery_data <- create_flat_recovery_data(true_params_df, recovered_params_df, "group")
    write.csv(recovery_data, recovery_csv_file, row.names = FALSE)
    
    # Verify file was created successfully
    if (!file.exists(recovery_csv_file)) {
      stop("Recovery CSV file was not created successfully")
    }
    
    message("Recovery data successfully saved to: ", recovery_csv_file)
  }, error = function(e) {
    warning("Failed to write recovery data: ", e$message)
  })
} else {
  # Hierarchical model analysis
  true_param_vals <- sub_params %>% 
    ungroup() %>% 
    as.data.frame() %>% 
    select(-c("idx"))
  
  recovered_params <- extract_recovery_parameters(
    fit = fits, # Could be a single fit object for hierarchical
    model_type = "hierarchical", 
    parameters = model_params,
    summary_fn = median
  )
  
  recovered_params = recovered_params$individual %>%
    pivot_wider(names_from = parameter, 
                values_from = value,
                id_cols = subject) %>%
    select(-subject)
  
  # Create CSV with recovery data with improved error handling - handle hierarchical structure
  tryCatch({
    # Validate that data exists
    if (ncol(true_param_vals) == 0 || nrow(true_param_vals) == 0) {
      stop("Empty true parameters data frame")
    }
    if (is.null(recovered_params)) {
      stop("Invalid recovered parameters structure")
    }
    
    # Ensure output directory exists
    output_dir_path <- dirname(recovery_csv_file)
    if (!dir.exists(output_dir_path)) {
      dir.create(output_dir_path, recursive = TRUE)
    }
    
    recovery_data <- create_flat_recovery_data(as.data.frame(true_param_vals), as.data.frame(recovered_params), "group")
    write.csv(recovery_data, recovery_csv_file, row.names = FALSE)
    
    # Verify file was created successfully
    if (!file.exists(recovery_csv_file)) {
      stop("Recovery CSV file was not created successfully")
    }
    
    message("Recovery data successfully saved to: ", recovery_csv_file)
  }, error = function(e) {
    warning("Failed to write recovery data: ", e$message)
  })
}

# Generate RMD analysis file with improved error handling
output_rmd <- file.path(
  output_rec_dir,
  generate_bids_filename(
    prefix = NULL,
    task = task,
    group = group_type,
    model = model_name,
    cohort = cohort,
    ses = session,
    additional_tags = list(
      "type" = "rec",
      "desc" = "analysis"
    ),
    ext = "Rmd"
  )
)

tryCatch({
  # Verify input file exists
  if (!file.exists(recovery_csv_file)) {
    stop("Recovery CSV file not found: ", recovery_csv_file)
  }
  
  # Create output directory if needed
  rmd_dir <- dirname(output_rmd)
  if (!dir.exists(rmd_dir)) {
    dir.create(rmd_dir, recursive = TRUE)
  }
  
  # Generate RMD file
  message("\nGenerating analysis RMD file...")
  generated_file <- generate_recovery_rmd(recovery_csv_file, output_rmd, opt$task, opt$model, opt$group, opt$render)
  
  # Verify the file was created
  if (!file.exists(generated_file)) {
    stop("RMD file generation failed: File not found")
  }
  
  message("RMD analysis file generated: ", generated_file)
}, error = function(e) {
  warning("Failed to generate RMD file: ", e$message, "\n",
          "You can generate it manually by running:\n",
          "Rscript analyze_recovery_results.R --input ", recovery_csv_file)
})

# Show output paths
message("\nParameter recovery data extraction complete.")
message("- Combined fit file saved to: ", combined_fit_file)
message("- Recovery data CSV saved to: ", recovery_csv_file)

# Explicitly exit with success status
quit(status = 0, save = "no")