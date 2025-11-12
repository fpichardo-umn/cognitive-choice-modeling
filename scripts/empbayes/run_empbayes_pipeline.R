#!/usr/bin/env Rscript

# Empirical Bayes Pipeline Orchestrator
# This script coordinates the execution of all empirical Bayes steps

suppressPackageStartupMessages({
  library(optparse)
  library(here)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--step"), type="character", default="all", 
              help="Which step(s) to run: 'all', '1', '2', '3', '4', or ranges like '1-3' (default: all)"),
  make_option(c("-k", "--task"), type="character", default=NULL, 
              help="Task name"),
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="Model name"),
  make_option(c("-s", "--source"), type="character", default=NULL, 
              help="Data source (cohort)"),
  make_option(c("--ses"), type="character", default=NULL, 
              help="Session identifier (optional)"),
  make_option(c("--model_status"), type="character", default=NULL,
              help="Model status (canonical/experimental/working)"),
  
  # Step 1: Subject selection
  make_option(c("--n_hier"), type="integer", default=NULL, 
              help="Number of subjects for hierarchical model"),
  make_option(c("--hier_subs_file"), type="character", default=NULL, 
              help="Pre-specified hierarchical subject list"),
  
  # Step 4: Empirical Bayes fitting
  make_option(c("--subjects"), type="character", default="all", 
              help="Subjects to fit (e.g., '1-100', 'all')"),
  make_option(c("--n_subs"), type="integer", default=NULL, 
              help="Number of subjects to fit"),
  
  # Data quality parameters (used in steps 1 and 4)
  make_option(c("--n_trials"), type="integer", default=120, 
              help="Minimum number of trials"),
  make_option(c("--RTbound_min_ms"), type="integer", default=50, 
              help="RT minimum bound in milliseconds"),
  make_option(c("--RTbound_max_ms"), type="integer", default=2500, 
              help="RT maximum bound in milliseconds"),
  make_option(c("--rt_method"), type="character", default="remove", 
              help="RT preprocessing method"),
  
  # Fitting parameters (used in steps 2 and 4)
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
              help="Checkpoint interval"),
  
  # Fitting method parameters
  make_option(c("--fitting_method"), type="character", default="mcmc",
              help="Fitting method for hierarchical model: mcmc or pathfinder (default: mcmc)"),
  make_option(c("--pf_num_paths"), type="integer", default=4,
              help="Pathfinder: number of paths (default: 4)"),
  make_option(c("--pf_draws"), type="integer", default=1000,
              help="Pathfinder: number of final draws (default: 1000)"),
  make_option(c("--pf_single_path_draws"), type="integer", default=250,
              help="Pathfinder: draws per path (default: 250)"),
  
  # Other parameters
  make_option(c("--seed"), type="integer", default=29518, 
              help="Random seed"),
  make_option(c("--parallel"), action="store_true", default=FALSE, 
              help="Use parallel processing in step 4"),
  make_option(c("--cores"), type="integer", default=4, 
              help="Number of cores for parallel processing"),
  make_option(c("--dry_run"), action="store_true", default=FALSE, 
              help="Perform a dry run")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Options used:\n")
dput(opt)

# Check required arguments
if (is.null(opt$task) || is.null(opt$model) || is.null(opt$source)) {
  stop("Task, model, and source are required. Use -k, -m, and -s options.")
}

# Parse step specification
parse_step_spec <- function(spec) {
  spec <- tolower(spec)
  
  if (spec == "all") {
    return(1:4)
  }
  
  # Handle ranges like "1-3"
  if (grepl("-", spec)) {
    parts <- as.integer(strsplit(spec, "-")[[1]])
    return(parts[1]:parts[2])
  }
  
  # Handle comma-separated like "1,3,4"
  if (grepl(",", spec)) {
    return(as.integer(strsplit(spec, ",")[[1]]))
  }
  
  # Single step
  return(as.integer(spec))
}

steps_to_run <- parse_step_spec(opt$step)

# Validate steps
if (any(steps_to_run < 1) || any(steps_to_run > 4)) {
  stop("Invalid step specification. Steps must be between 1 and 4.")
}

# Get script directory
script_dir <- file.path(here::here(), "scripts", "empbayes")

# Function to run a script with arguments
run_script <- function(script_name, args) {
  script_path <- file.path(script_dir, script_name)
  
  if (!file.exists(script_path)) {
    stop("Script not found: ", script_path)
  }
  
  cmd <- c(script_path, args)
  cat("\n", strrep("=", 80), "\n")
  cat("Executing:", script_name, "\n")
  cat(strrep("=", 80), "\n")
  cat("Command: Rscript", paste(cmd, collapse=" "), "\n\n")
  
  result <- system2("Rscript", cmd)
  
  if (result != 0) {
    stop("Script failed with exit code: ", result)
  }
  
  cat("\n", strrep("=", 80), "\n")
  cat("Completed:", script_name, "\n")
  cat(strrep("=", 80), "\n\n")
  
  return(result)
}

# Build common arguments
common_args <- c(
  "-k", opt$task,
  "-m", opt$model,
  "-s", opt$source
)

if (!is.null(opt$ses)) {
  common_args <- c(common_args, "--ses", opt$ses)
}

if (!is.null(opt$model_status)) {
  common_args <- c(common_args, "--model_status", opt$model_status)
}

common_args <- c(common_args, "--seed", as.character(opt$seed))

if (opt$dry_run) {
  common_args <- c(common_args, "--dry_run")
}

# Print pipeline configuration
cat("\n")
cat(strrep("=", 80), "\n")
cat("EMPIRICAL BAYES PIPELINE\n")
cat(strrep("=", 80), "\n")
cat("Task:", opt$task, "\n")
cat("Model:", opt$model, "\n")
cat("Source:", opt$source, "\n")
if (!is.null(opt$ses)) cat("Session:", opt$ses, "\n")
cat("Steps to run:", paste(steps_to_run, collapse=", "), "\n")
cat("Dry run:", opt$dry_run, "\n")
cat(strrep("=", 80), "\n\n")

# Step 1: Select subjects for hierarchical model
if (1 %in% steps_to_run) {
  step1_args <- c(
    common_args,
    "--n_trials", as.character(opt$n_trials),
    "--RTbound_min_ms", as.character(opt$RTbound_min_ms),
    "--RTbound_max_ms", as.character(opt$RTbound_max_ms),
    "--rt_method", opt$rt_method
  )
  
  if (!is.null(opt$n_hier)) {
    step1_args <- c(step1_args, "--n_hier", as.character(opt$n_hier))
  }
  
  if (!is.null(opt$hier_subs_file)) {
    step1_args <- c(step1_args, "--hier_subs_file", opt$hier_subs_file)
  }
  
  run_script("select_empbayes_subjects.R", step1_args)
}

# Step 2: Fit hierarchical model
if (2 %in% steps_to_run) {
  step2_args <- c(
    common_args,
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
    "--fitting_method", opt$fitting_method,
    "--pf_num_paths", as.character(opt$pf_num_paths),
    "--pf_draws", as.character(opt$pf_draws),
    "--pf_single_path_draws", as.character(opt$pf_single_path_draws)
  )
  
  if (!is.null(opt$hier_subs_file)) {
    step2_args <- c(step2_args, "--hier_subs_file", opt$hier_subs_file)
  }
  
  run_script("fit_empbayes_hierarchical.R", step2_args)
}

# Step 3: Generate informative priors
if (3 %in% steps_to_run) {
  step3_args <- common_args
  
  run_script("generate_empbayes_priors.R", step3_args)
}

# Step 4: Fit individual subjects with empirical Bayes priors
if (4 %in% steps_to_run) {
  step4_args <- c(
    common_args,
    "--n_trials", as.character(opt$n_trials),
    "--RTbound_min_ms", as.character(opt$RTbound_min_ms),
    "--RTbound_max_ms", as.character(opt$RTbound_max_ms),
    "--rt_method", opt$rt_method,
    "--n_warmup", as.character(opt$n_warmup),
    "--n_iter", as.character(opt$n_iter),
    "--n_chains", as.character(opt$n_chains),
    "--adapt_delta", as.character(opt$adapt_delta),
    "--max_treedepth", as.character(opt$max_treedepth),
    "--check_iter", as.character(opt$check_iter)
  )
  
  if (!is.null(opt$subjects)) {
    step4_args <- c(step4_args, "--subjects", opt$subjects)
  }
  
  if (!is.null(opt$n_subs)) {
    step4_args <- c(step4_args, "--n_subs", as.character(opt$n_subs))
  }
  
  if (opt$parallel) {
    step4_args <- c(step4_args, "--parallel", "--cores", as.character(opt$cores))
  }
  
  run_script("fit_empbayes_batch.R", step4_args)
}

# Pipeline complete
cat("\n")
cat(strrep("=", 80), "\n")
cat("PIPELINE COMPLETE\n")
cat(strrep("=", 80), "\n")
cat("All requested steps completed successfully.\n")
cat(strrep("=", 80), "\n\n")
