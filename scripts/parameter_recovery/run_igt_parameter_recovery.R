#!/usr/bin/env Rscript

# Parameter recovery test script for standard IGT
# This script simulates data from the standard IGT and tests parameter recovery

library(R6)
library(data.table)
library(yaml)
library(optparse)
library(tidyverse)
library(here)

# Source necessary code files
base_dir <- here::here("scripts", "parameter_recovery")
source(file.path(base_dir, "helper_functions_PR.R"))

# Simulation code
sim_dir <- here::here("scripts", "simulation")
source(file.path(sim_dir, "tasks", "base_task.R"))
source(file.path(sim_dir, "tasks", "igt", "igt_task.R"))
source(file.path(sim_dir, "models", "base_model.R"))
source(file.path(sim_dir, "models", "igt", "igt_ev_model.R"))
source(file.path(sim_dir, "simulator.R"))
source(file.path(sim_dir, "param_gen.R"))

# Define command-line options
option_list <- list(
  make_option(c("-m", "--model"), type = "character", default = "ev",
              help = "Model to use for parameter recovery [default: %default]"),
  make_option(c("-k", "--task"), type = "character", default = "igt",
              help = "Task to simulate [default: %default]"),
  make_option(c("-n", "--n_subjects"), type = "integer", default = 100,
              help = "Number of subjects to simulate [default: %default]"),
  make_option(c("-t", "--n_trials"), type = "integer", default = 100,
              help = "Number of trials per subject [default: %default]"),
  make_option(c("-b", "--n_blocks"), type = "integer", default = 5,
              help = "Number of blocks [default: %default]"),
  make_option(c("--trials_per_block"), type = "integer", default = 20,
              help = "Trials per block [default: %default]"),
  make_option(c("-o", "--output_dir"), type = "character", default = "results",
              help = "Output directory [default: %default]"),
  make_option(c("-s", "--seed"), type = "integer", default = 42,
              help = "Random seed [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Create task instance
task <- IGTTask$new()

# Load model configuration
model_config_path <- file.path(sim_dir, "config", "tasks", opts$task, "models", paste0(opts$model, ".yaml"))
model_config <- yaml::read_yaml(model_config_path)

# Create model instance
model <- IGTEVModel$new(task)

# Generate parameters
param_info <- model$get_parameter_info()
n_subjects <- opts$n_subjects

# Generate random parameters for simulation
params <- generate_random_parameters(
  parameter_info = param_info,
  n_subjects = n_subjects,
  seed = opts$seed
)

# Add simulation index
params$idx <- "sim_params"

# Simulate data
sim_data <- simulate_data(
  task = task,
  model = model,
  parameters = params,
  n_trials = opts$n_trials,
  n_blocks = opts$n_blocks,
  trials_per_block = opts$trials_per_block,
  seed = opts$seed,
  return_format = "data.table"
)

# Print summary of simulated data
cat("Simulated", nrow(sim_data), "trials for", n_subjects, "subjects\n")
cat("Mean outcome:", mean(sim_data$outcome, na.rm = TRUE), "\n")

# Create deck selection frequency summary
deck_freq <- sim_data %>%
  group_by(param_idx) %>%
  summarize(
    A_freq = mean(choice == 1, na.rm = TRUE),
    B_freq = mean(choice == 2, na.rm = TRUE),
    C_freq = mean(choice == 3, na.rm = TRUE),
    D_freq = mean(choice == 4, na.rm = TRUE)
  )

cat("\nDeck selection frequencies:\n")
print(colMeans(deck_freq[, -1]))

# Save results
output_dir <- file.path(here::here(), "Data", "sim", "igt_recovery")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(
  output_dir,
  paste0(
    "igt_sim_data_",
    opts$model,
    "_n", n_subjects,
    "_t", opts$n_trials,
    ".rds"
  )
)

saveRDS(
  list(
    data = sim_data,
    parameters = params,
    config = list(
      model = opts$model,
      task = opts$task,
      n_subjects = n_subjects,
      n_trials = opts$n_trials,
      n_blocks = opts$n_blocks,
      trials_per_block = opts$trials_per_block,
      seed = opts$seed
    )
  ),
  file = output_file
)

cat("\nResults saved to:", output_file, "\n")