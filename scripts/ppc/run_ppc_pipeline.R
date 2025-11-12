#!/usr/bin/env Rscript

#' Run Complete PPC Pipeline
#' @description Execute all steps of the PPC analysis pipeline

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(dplyr)
})

# ---- Define command line options ----
option_list = list(
  make_option(c("-m", "--model"), type="character", help="Model name (e.g., ev, pvldelta)"),
  make_option(c("-k", "--task"), type="character", default="igt_mod", help="Task name"),
  make_option(c("-c", "--cohort"), type="character", help="Cohort identifier"),
  make_option(c("--ses"), type="character", default=NULL, help="Session identifier (optional)"),
  make_option(c("-g", "--group"), type="character", default="sing", help="Fit type: sing (individual) or hier (hierarchical)"),
  make_option(c("--group_name"), type="character", default="batch_001", help="Batch identifier (only used for individual fits)"),
  make_option(c("-f", "--fit_file"), type="character", default=NULL, 
              help="Path to fit file (optional, built from task, group, model if not provided)"),
  make_option(c("-n", "--n_sims"), type="integer", default=100, 
              help="Number of simulations per subject"),
  make_option(c("-b", "--block_size"), type="integer", default=20, 
              help="Number of trials per block"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory (optional)"),
  make_option(c("-p", "--parallel"), action="store_true", default=FALSE, 
              help="Use parallel processing"),
  make_option(c("--n_cores"), type="integer", default=2, 
              help="Number of cores for parallel processing"),
  make_option(c("-s", "--steps"), type="character", default="all", 
              help="Steps to run (comma-separated list): simulate,stats,loglik,report,all"),
  make_option(c("--force"), action="store_true", default=FALSE, 
              help="Force re-running of steps even if output exists"),
  make_option(c("--render"), action="store_true", default=FALSE,
              help="Render the Rmd file to HTML after generating it"),
  make_option(c("--sampling"), type="character", default="weighted", 
              help="Sampling method for posterior draws: random, width, or weighted"),
  make_option(c("--width_control"), type="numeric", default=0.95, 
              help="Width control parameter for width sampling method (0-1)"),
  make_option(c("--rt_method"), type="character", default="remove", 
              help="RT handling method: all, remove, force, or adaptive"),
  make_option(c("--RTbound_min_ms"), type="numeric", default=100, 
              help="RT lower bound in milliseconds"),
  make_option(c("--RTbound_max_ms"), type="numeric", default=2500, 
              help="RT upper bound in milliseconds"),
  make_option(c("--exclude_file"), type="character", default=NULL, 
              help="Path to file with subject IDs to exclude"),
  make_option(c("--ic_method"), type="character", default="loo", 
              help="Information criterion method: loo or waic"),
  make_option(c("--n_samples_loglik"), type="character", default="all",
              help="Number of posterior samples for loglik (default: all, or integer)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$model) || is.null(opt$cohort)) {
  stop("Model name and cohort are required. Use --model and --cohort flags.")
}

# Validate group parameter
if (!opt$group %in% c("sing", "hier")) {
  stop("Invalid --group parameter. Must be 'sing' (individual) or 'hier' (hierarchical).")
}

# Set up environment - import helper modules
script_dir <- file.path(here::here(), "scripts")
source(file.path(script_dir, "ppc", "helpers", "helper_ppc_dirs.R"))
source(file.path(script_dir, "helpers", "helper_functions_cmdSR.R"))

# Ensure PPC directories exist for the given task
directory_info <- ensure_ppc_dirs(opt$task, opt$cohort, opt$ses)

# Set output directory
output_dir <- if(!is.null(opt$output_dir)) opt$output_dir else get_ppc_dir(opt$task, opt$cohort, opt$ses)

# Determine which steps to run
steps_to_run <- opt$steps
if (steps_to_run == "all") {
  # Run loglik before stats since:
  # 1. Both only depend on simulation (independent)
  # 2. Loglik is fast, stats is slow
  # 3. Get IC results quickly while stats runs
  steps_to_run <- c("simulate", "loglik", "stats", "report")
} else {
  steps_to_run <- strsplit(steps_to_run, ",")[[1]]
}

# Check for valid steps
valid_steps <- c("simulate", "stats", "loglik", "report", "all")
invalid_steps <- setdiff(steps_to_run, valid_steps)
if (length(invalid_steps) > 0) {
  stop("Invalid steps specified: ", paste(invalid_steps, collapse = ", "), 
       ". Valid options are: ", paste(valid_steps, collapse = ", "))
}

# Get file paths for outputs using helper functions
# Output file naming should reflect the source of parameters:
# - For hierarchical fits: use "hier" (even though we simulate individuals)
# - For individual fits: use group_name (batch identifier)
group_for_output_files <- if(opt$group == "hier") "hier" else opt$group_name
sim_file <- get_ppc_sim_file_path(opt$task, opt$model, group_for_output_files, opt$cohort, opt$ses)
stats_file <- file.path(get_ppc_stats_dir(opt$task, opt$cohort, opt$ses),
                        generate_bids_filename("ppc_summary", opt$task, group_for_output_files, opt$model, 
                                               "csv", cohort = opt$cohort, ses = opt$ses))
loglik_file <- get_ppc_loglik_file_path(opt$task, opt$model, group_for_output_files, opt$cohort, opt$ses)
report_file <- get_ppc_report_file_path(opt$task, opt$model, group_for_output_files, opt$cohort, opt$ses)

# Create command line arguments for each script
run_simulation_args <- paste0(
  " --model ", opt$model,
  " --task ", opt$task,
  " --cohort ", opt$cohort,
  " --group ", opt$group,
  " --group_name ", opt$group_name,
  " --n_sims ", opt$n_sims,
  " --sampling ", opt$sampling,
  " --width_control ", opt$width_control,
  " --rt_method ", opt$rt_method,
  " --RTbound_min_ms ", opt$RTbound_min_ms,
  " --RTbound_max_ms ", opt$RTbound_max_ms
)

# Add session if provided
if (!is.null(opt$ses)) {
  run_simulation_args <- paste0(run_simulation_args, " --ses ", opt$ses)
}

if (!is.null(opt$fit_file)) {
  run_simulation_args <- paste0(run_simulation_args, " --fit_file \"", opt$fit_file, "\"")
}

if (opt$parallel) {
  run_simulation_args <- paste0(run_simulation_args, " --parallel --n_cores ", opt$n_cores)
}

run_stats_args <- paste0(
  " --model ", opt$model,
  " --task ", opt$task,
  " --cohort ", opt$cohort,
  " --group ", group_for_output_files,  # Use group_for_output_files for output file naming
  " --block_size ", opt$block_size,
  " --sim_file \"", sim_file, "\""
)

# Add session if provided
if (!is.null(opt$ses)) {
  run_stats_args <- paste0(run_stats_args, " --ses ", opt$ses)
}

if (!is.null(opt$exclude_file)) {
  run_simulation_args <- paste0(run_simulation_args, " --exclude_file \"", opt$exclude_file, "\"")
  run_stats_args <- paste0(run_stats_args, " --exclude_file \"", opt$exclude_file, "\"")
}

# Build loglik args - ALWAYS use fit_file (LOOIC should use all posterior draws)
# Determine fit file path
if (!is.null(opt$fit_file)) {
  fit_file_to_use <- opt$fit_file
} else {
  # Construct fit file path
  fit_dir <- get_fits_output_dir(opt$task, "fit", opt$cohort, opt$ses)
  group_identifier <- if(opt$group == "hier") "hier" else opt$group_name
  fit_file_to_use <- file.path(fit_dir, 
                               generate_bids_filename(NULL, opt$task, group_identifier, opt$model,
                                                     ext = "rds", cohort = opt$cohort, ses = opt$ses,
                                                     additional_tags = list("type" = "fit", "desc" = "output")))
}

# Build loglik args with user-specified n_samples
run_loglik_args <- paste0(
  " --model ", opt$model,
  " --task ", opt$task,
  " --cohort ", opt$cohort,
  " --group ", group_for_output_files,
  " --fit_file \"", fit_file_to_use, "\"",
  " --ic_method ", opt$ic_method,
  " --rt_method ", opt$rt_method,
  " --RTbound_min_ms ", opt$RTbound_min_ms,
  " --RTbound_max_ms ", opt$RTbound_max_ms,
  " --n_samples ", opt$n_samples_loglik
)
message("Log-likelihood will use ", opt$n_samples_loglik, " posterior draws from fit file")

# Add exclude_file to loglik args if provided
if (!is.null(opt$exclude_file)) {
  run_loglik_args <- paste0(run_loglik_args, " --exclude_file \"", opt$exclude_file, "\"")
}

# Add session if provided
if (!is.null(opt$ses)) {
  run_loglik_args <- paste0(run_loglik_args, " --ses ", opt$ses)
}

run_report_args <- paste0(
  " --model ", opt$model,
  " --task ", opt$task,
  " --cohort ", opt$cohort,
  " --group ", group_for_output_files,  # Use group_for_output_files for output file naming
  " --stats_file \"", stats_file, "\""
)

# Add session if provided
if (!is.null(opt$ses)) {
  run_report_args <- paste0(run_report_args, " --ses ", opt$ses)
}

if (opt$force) {
  run_report_args <- paste0(run_report_args, " --force")
}

if (opt$render) {
  run_report_args <- paste0(run_report_args, " --render")
}

# ---- Step 1: Simulation ----
if ("simulate" %in% steps_to_run) {
  if (file.exists(sim_file) && !opt$force) {
    message("Simulation output already exists. Skipping simulation step. Use --force to override.")
  } else {
    message("Step 1: Running simulations")
    simulation_script <- file.path(script_dir, "ppc", "run_simulation.R")
    cmd <- paste0("Rscript \"", simulation_script, "\"", run_simulation_args)
    message("Executing command: ", cmd)
    system(cmd)
    
    if (!file.exists(sim_file)) {
      stop("Simulation failed. Output file not found: ", sim_file)
    }
  }
} else {
  message("Skipping simulation step as requested.")
}

# ---- Step 2: Calculate log-likelihood (independent of simulation!) ----
if ("loglik" %in% steps_to_run) {
  if (file.exists(loglik_file) && !opt$force) {
    message("Log-likelihood output already exists. Skipping log-likelihood step. Use --force to override.")
  } else {
    # Loglik can run with or without simulation file
    # It will use sim_file if available, otherwise load data directly from fit_file
    
    message("Step 2: Calculating log-likelihood")
    loglik_script <- file.path(script_dir, "ppc", "run_loglik.R")
    cmd <- paste0("Rscript \"", loglik_script, "\"", run_loglik_args)
    message("Executing command: ", cmd)
    system(cmd)
    
    if (!file.exists(loglik_file)) {
      stop("Log-likelihood calculation failed. Output file not found: ", loglik_file)
    }
  }
} else {
  message("Skipping log-likelihood step as requested.")
}

# ---- Step 3: Calculate statistics (can be slow, so runs after loglik) ----
if ("stats" %in% steps_to_run) {
  if (file.exists(stats_file) && !opt$force) {
    message("Statistics output already exists. Skipping statistics step. Use --force to override.")
  } else {
    # Check if simulation file exists
    if (!file.exists(sim_file)) {
      stop("Simulation output not found. Run simulation step first: ", sim_file)
    }
    
    message("Step 3: Calculating statistics")
    stats_script <- file.path(script_dir, "ppc", "run_stats.R")
    cmd <- paste0("Rscript \"", stats_script, "\"", run_stats_args)
    message("Executing command: ", cmd)
    system(cmd)
    
    if (!file.exists(stats_file)) {
      stop("Statistics calculation failed. Output file not found: ", stats_file)
    }
  }
} else {
  message("Skipping statistics step as requested.")
}



# ---- Step 4: Generate report ----
if ("report" %in% steps_to_run) {
  if (file.exists(report_file) && !opt$force) {
    message("Report already exists. Skipping report generation. Use --force to override.")
  } else {
    # Check if statistics file exists
    if (!file.exists(stats_file)) {
      stop("Statistics output not found. Run statistics step first: ", stats_file)
    }
    
    message("Step 4: Generating report")
    report_script <- file.path(script_dir, "ppc", "generate_ppc_report.R")
    cmd <- paste0("Rscript \"", report_script, "\"", run_report_args)
    message("Executing command: ", cmd)
    system(cmd)
    
    if (!file.exists(report_file)) {
      warning("Report generation may have failed. Output file not found: ", report_file)
    }
  }
} else {
  message("Skipping report generation as requested.")
}

message("PPC pipeline execution complete.")
message("Output files:")
message("  Simulations: ", sim_file)
message("  Statistics: ", stats_file)
message("  Log-likelihood: ", loglik_file)
message("  Report: ", report_file)