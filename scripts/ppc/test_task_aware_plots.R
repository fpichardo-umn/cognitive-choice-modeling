#!/usr/bin/env Rscript

# Test script for task-aware PPC visualizations
# This script tests the new task-aware functionality

# Load required libraries
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
})

# Source the updated functions
source(file.path(here::here(), "scripts", "ppc", "helpers", "task_config.R"))
source(file.path(here::here(), "scripts", "ppc", "visualization_functions.R"))

# Test 1: Get task configurations
cat("=== Testing Task Configurations ===\n")

# Test IGT configuration
igt_config <- get_task_config("igt")
cat("IGT task type:", igt_config$type, "\n")
cat("IGT applicable stats:", length(unlist(igt_config$statistics)), "total\n")

# Test mIGT configuration
migt_config <- get_task_config("igt_mod")
cat("mIGT task type:", migt_config$type, "\n")
cat("mIGT applicable stats:", length(unlist(migt_config$statistics)), "total\n")

# Test 2: Test available stat groups function
cat("\n=== Testing Available Stat Groups ===\n")

# Mock some available statistics for IGT
igt_available_stats <- c("deck1_freq", "deck2_freq", "good_deck_freq", "net_score", "win_stay", "rt_mean")
igt_groups <- get_available_stat_groups("igt", igt_available_stats)
cat("IGT available groups:", paste(igt_groups, collapse = ", "), "\n")

# Mock some available statistics for mIGT
migt_available_stats <- c("play_ratio", "play_ratio_deck1", "good_play_ratio", "net_score", "rt_mean_play")
migt_groups <- get_available_stat_groups("igt_mod", migt_available_stats)
cat("mIGT available groups:", paste(migt_groups, collapse = ", "), "\n")

# Test 3: Test validation function
cat("\n=== Testing Statistic Validation ===\n")

# Test some statistics
test_stats <- c("play_ratio", "deck1_freq", "win_stay", "rt_mean_play")
for (stat in test_stats) {
  igt_valid <- validate_statistic_for_task(stat, "igt")
  migt_valid <- validate_statistic_for_task(stat, "igt_mod")
  cat(sprintf("%-15s: IGT=%s, mIGT=%s\n", stat, igt_valid, migt_valid))
}

cat("\n=== Task-Aware PPC Functions Test Complete ===\n")
cat("If no errors appeared above, the basic