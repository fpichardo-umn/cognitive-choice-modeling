#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(optparse)
  library(here)
})

# Source helper functions
source(file.path(here::here(), "scripts", "parameter_recovery", "helper_functions_PR.R"))
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))

# Define command line options (since we're keeping these in each script)
option_list = list(
  make_option(c("-m", "--model"), type="character", help="Model name (e.g., ev, pvl)"),
  make_option(c("-t", "--task"), type="character", help="Task name (e.g., igt_mod)"),
  make_option(c("-g", "--group"), type="character", default="sing", help="Group type"),
  make_option(c("--cohort"), type="character", default=NULL,
              help="Cohort identifier (e.g., ahrb, es)"),
  make_option(c("--session"), type="character", default=NULL,
              help="Session identifier [default: %default]"),
  make_option(c("-p", "--param_file"), type="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="./Data/sim/txt/"),
  make_option(c("-b", "--n_blocks"), type="integer", default=6),
  make_option(c("-k", "--trials_per_block"), type="integer", default=20),
  make_option(c("-s", "--seed"), type="integer", default=12345),
  make_option(c("--RTbound_min_ms"), type="numeric", default=50, 
              help="RT lower bound in milliseconds [default: %default]"),
  make_option(c("--RTbound_max_ms"), type="numeric", default=1200000, 
              help="RT upper bound in milliseconds [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Get directory structure
dirs <- setup_directories(opt$task)

# Source required files
source_required_files(dirs$SIM_DIR, opt$task)
source(file.path(dirs$SIM_DIR, "simulator.R"))

# Extract model name and create full model name
model_name <- opt$model
task_name <- opt$task
group_type <- opt$group
cohort <- opt$cohort
session <- opt$session
full_model_name <- paste(task_name, group_type, model_name, sep="_")

# Load parameters
if (is.null(opt$param_file)){
  filename <- file.path(
    get_simulation_output_dir(opt$task, "parameters"),
    generate_bids_filename(
      prefix = NULL,
      task = opt$task,
      group = opt$group,
      model = opt$model,
      cohort = opt$cohort,
      ses = opt$session,
      additional_tags = list(
        "type" = "params",
        "desc" = "*",
        "n" = "*"
      ),
      ext = "rds"
    )
  )
  param_file = Sys.glob(filename)[1]
} else {
  param_file = opt$param_file
}
params <- readRDS(param_file)

# Initialize task and model
task <- initialize_task(task_name, dirs$SIM_DIR)

if (grepl("batch", opt$group)){
  group_type = "sing"
} else {
  group_type = opt$group
}
model <- initialize_model(model_name, task_name, task, dirs$SIM_DIR, group_type)

# Build task_params from command line RT bounds
task_params <- list(
  RTbound_min = opt$RTbound_min_ms / 1000,  # Convert ms to seconds
  RTbound_max = opt$RTbound_max_ms / 1000
)

cat("Using RT bounds: [", task_params$RTbound_min, ", ", 
    task_params$RTbound_max, "] seconds\n")

# Run simulation
sim_data <- simulate_data(
  task = task,
  model = model,
  parameters = params,
  n_trials = opt$n_blocks * opt$trials_per_block,
  n_blocks = opt$n_blocks,
  trials_per_block = opt$trials_per_block,
  seed = opt$seed,
  task_params = task_params
)

# Save simulated data as CSV
csv_filename <- file.path(
  get_simulation_output_dir(opt$task, "data"),
  "txt",
  generate_bids_filename(
    prefix = NULL,
    task = task_name,
    group = group_type,
    model = model_name,
    cohort = cohort,
    ses = session,
    additional_tags = list(
      "type" = "sim",
      "desc" = "data"
    ),
    ext = "csv"
  )
)
write.csv(data.frame(sim_data), csv_filename, row.names = FALSE)

# Save simulated data as RDS
rds_filename <- file.path(
  get_simulation_output_dir(opt$task, "data"),
  "rds",
  generate_bids_filename(
    prefix = NULL,
    task = task_name,
    group = group_type,
    model = model_name,
    cohort = cohort,
    ses = session,
    additional_tags = list(
      "type" = "sim",
      "desc" = "data"
    ),
    ext = "rds"
  )
)
saveRDS(sim_data, rds_filename)

cat("Simulation data saved to:", "\n")
cat(rds_filename)
