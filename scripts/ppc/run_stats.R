#!/usr/bin/env Rscript

#' Run Statistics for Posterior Predictive Checks
#' @description Calculate statistics and PPP values from observed and simulated data

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(foreign)  # For read.spss
})

# Define command line options
option_list = list(
  make_option(c("-m", "--model"), type="character", default=NULL, help="Model name"),
  make_option(c("-k", "--task"), type="character", default=NULL, help="Task name"),
  make_option(c("-c", "--cohort"), type="character", default=NULL, help="Cohort identifier"),
  make_option(c("-g", "--group"), type="character", default="batch_001", help="Group identifier for file naming"),
  make_option(c("--ses"), type="character", default=NULL, help="Session identifier"),
  make_option(c("-b", "--block_size"), type="integer", default=20, 
              help="Number of trials per block"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory (optional)"),
  make_option(c("-e", "--exclude_file"), type="character", default=NULL, 
              help="File with subject IDs to exclude"),
  make_option(c("-s", "--sim_file"), type="character", default=NULL, 
              help="Path to simulation results file (optional)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Options used:\n")
dput(opt)

# Check required arguments
if (is.null(opt$model) || is.null(opt$task) || is.null(opt$cohort)) {
  stop("Model name, task name, and cohort are all required.")
}

# Source helper functions - use unified helper entry point
script_dir <- file.path(here::here(), "scripts")
source(file.path(script_dir, "ppc", "helpers", "helper_functions_ppc.R"))
source(file.path(script_dir, "ppc", "statistics_functions.R"))  # For backward compatibility

# Ensure PPC directories exist
directory_info <- ensure_ppc_dirs(opt$task, opt$cohort, opt$ses)

# Set output directory
output_dir <- if(!is.null(opt$output_dir)) opt$output_dir else get_ppc_stats_dir(opt$task, opt$cohort, opt$ses)

# Determine model type based on model name
model_type <- if(opt$model == "ddm") {
  "SSM"  # Pure Sequential Sampling Model
} else if(grepl("ddm", opt$model)) {
  "RL_SSM"  # Combined RL and Sequential Sampling Model
} else {
  "RL"  # Pure Reinforcement Learning Model
}
message("Determined model type: ", model_type)

# Build or use simulation file path
if (is.null(opt$sim_file)) {
  # Build standard simulation file path using helper function
  sim_file <- get_ppc_sim_file_path(opt$task, opt$model, opt$group, opt$cohort, opt$ses)
  message("Using standard simulation file path: ", sim_file)
} else {
  sim_file <- opt$sim_file
  message("Using provided simulation file: ", sim_file)
}

# Load simulation results with error handling
message("Loading simulation results...")
simulation_results <- load_rds_safe(sim_file, "PPC simulation results")
if (is.null(simulation_results)) {
  stop("Failed to load simulation results from: ", sim_file)
}
subject_ids <- names(simulation_results)
message("Found ", length(subject_ids), " subjects in simulation results")

# Load observed statistics
message("Loading observed statistics...")
observed_stats <- list()

# Try to load pre-computed observed statistics with dynamic paths
SAFE_DATA_DIR <- get_safe_data_dir()
# Build dynamic file paths based on task/cohort/session
session_stats_file <- file.path(SAFE_DATA_DIR, opt$cohort,
                                generate_bids_filename("subs_session_stats", 
                                                       opt$task, group = NULL, model = NULL, 
                                                       ext = "csv", cohort = opt$cohort, ses = opt$ses)
                                )
block_stats_file <- file.path(SAFE_DATA_DIR, opt$cohort,
                              generate_bids_filename("subs_block_stats", 
                                                     opt$task, group = NULL, model = NULL, 
                                                     ext = "csv", cohort = opt$cohort, ses = opt$ses)
                              )

# Check if pre-computed stats exist and can be loaded
stats_exist <- file.exists(session_stats_file) && file.exists(block_stats_file)
stats_loaded <- FALSE

if (stats_exist) {
  message("Loading pre-computed observed statistics...")
  observed_stats$session <- load_csv_safe(session_stats_file, "session statistics")
  observed_stats$blocks <- load_csv_safe(block_stats_file, "block statistics")
  
  if (!is.null(observed_stats$session) && !is.null(observed_stats$blocks)) {
    stats_loaded <- TRUE
  } else {
    message("Failed to load pre-computed stats, will regenerate from raw data")
  }
}

if (!stats_loaded) {
  # Fall back to generating from raw data using proper load_data function
  message("Pre-computed stats not found, generating from raw data...")
  
  # Load data using the same approach as simulation step
  source(file.path(script_dir, "helpers", "helper_load.R"))
  all_data <- load_data(opt$task, opt$cohort, opt$ses)
  
  # Filter observed data for subjects that have simulation results
  observed_data <- all_data[all_data$subjID %in% subject_ids, ]
  message("Using observed data for ", length(unique(observed_data$subjID)), " subjects")
  
  # Calculate observed statistics
  observed_stats_temp <- calculate_observed_statistics(observed_data, model_type, opt$block_size, opt$task, opt$cohort)
  
  # Align with expected format
  observed_stats <- create_observed_stats_from_temp(observed_stats_temp)
  
  # Save observed data with error handling
  save_csv_safe(observed_stats$session, session_stats_file, "session statistics")
  save_csv_safe(observed_stats$blocks, block_stats_file, "block statistics")
}

# Load exclude list if provided
exclude_subjects <- NULL
if (!is.null(opt$exclude_file) && file.exists(opt$exclude_file)) {
  exclude_subjects <- readLines(opt$exclude_file)
  message(paste("Excluding", length(exclude_subjects), "subjects"))
}

# Filter excluded subjects
observed_stats$session <- observed_stats$session %>%
  filter(!sid %in% exclude_subjects)

observed_stats$blocks <- observed_stats$blocks %>%
  filter(!sid %in% exclude_subjects)

simulation_results = simulation_results[setdiff(names(simulation_results), exclude_subjects)]

# Calculate simulation statistics
message("Calculating simulation statistics...")
simulation_stats <- calculate_simulation_statistics(
  simulation_results, model_type, opt$block_size, opt$task
)

# Calculate PPP statistics
message("Calculating posterior predictive p-values...")
ppc_stats <- calculate_ppc_statistics(observed_stats, simulation_stats)

# Save results with error handling
stats_file <- get_ppc_stats_file_path(opt$task, opt$model, opt$group, opt$cohort, opt$ses)

success <- save_rds_safe(ppc_stats, stats_file, "PPC statistics")
if (!success) {
  warning("Failed to save PPC statistics, but continuing...")
}

# Generate summary statistics
ppc_summary <- summarize_ppc_statistics(ppc_stats)

# Add task information
ppc_summary$task <- opt$task

# Save summary using helper function
summary_file <- save_ppc_statistics(ppc_summary, opt$task, opt$group, opt$model, opt$cohort, get_ppc_stats_dir(opt$task, opt$cohort, opt$ses), session = opt$ses)
message("Summary statistics saved to: ", summary_file)

# Basic report with extreme PPP values
extreme_summary <- ppc_summary %>%
  filter(extreme_ppp) %>%
  arrange(desc(abs(ppp - 0.5)))

extreme_count <- nrow(extreme_summary)
total_count <- nrow(ppc_summary)
extreme_percent <- 100 * extreme_count / total_count

message("\nPPC Analysis Summary:")
message(paste0(extreme_count, " of ", total_count, 
               " statistics (", round(extreme_percent, 1), 
               "%) show extreme PPP values (< 0.05 or > 0.95)"))

message("\nPPC analysis complete.")
