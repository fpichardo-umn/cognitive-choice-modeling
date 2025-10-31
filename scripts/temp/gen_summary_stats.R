#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(cmdstanr)
  library(posterior)
  library(foreign)
  library(dplyr)
  library(tidyr)
})

# Set up directories
PROJ_DIR <- here::here()
DATA_DIR <- file.path(PROJ_DIR, "Data")
SAFE_DATA_DIR <- file.path(DATA_DIR, "AHRB")
SCRIPT_DIR <- file.path(PROJ_DIR, "scripts")

# Load helper functions
helper_functions_path <- file.path(SCRIPT_DIR, "helper_functions_cmdSR.R")
if (!file.exists(helper_functions_path)) {
  stop("helper_functions.R not found. Expected path: ", helper_functions_path)
}
source(helper_functions_path)

# Check for data file existence
wave1.sav.file <- file.path(SAFE_DATA_DIR, "modigt_data_Wave1.sav")
if (!file.exists(wave1.sav.file)) {
  stop("Data file not found. Expected path: ", wave1.sav.file)
}

# Load data
wave1.raw <- read.spss(wave1.sav.file, to.data.frame = TRUE)

# Calculate metrics
choice_stats <- calculate_choice_metrics(wave1.raw)
rt_stats <- calculate_rt_metrics(wave1.raw)

# Combine statistics
participant_stats <- choice_stats %>%
  left_join(rt_stats, by = "sid")

# Create filename with wave number
wave_number <- as.numeric(gsub(".*Wave([0-9]+)\\.sav$", "\\1", wave1.sav.file))
filename <- paste0("participant_stats_wave", wave_number, ".csv")
filepath <- file.path(SAFE_DATA_DIR, filename)

# Save to CSV
write.csv(participant_stats, filepath, row.names = FALSE)

cat("Participant statistics saved to:", filepath, "\n")
