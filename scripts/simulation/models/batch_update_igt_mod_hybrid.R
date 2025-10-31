#!/usr/bin/env Rscript

# Batch update script for IGT_MOD Hybrid models
# Updates calculate_loglik to new format: (data, parameters, task_params)

library(stringr)

# List of files to update (remaining IGT_MOD hybrid models)
files_to_update <- c(
  # experimental models
  "igt_mod/experimental/igt_mod_ev_tic_ddm_model.R",
  "igt_mod/experimental/igt_mod_ev_tic_model.R",
  "igt_mod/experimental/igt_mod_pvldecay_ddm_b1b_model.R",
  "igt_mod/experimental/igt_mod_pvldecay_ddm_b1p2_model.R",
  "igt_mod/experimental/igt_mod_pvldecay_tic_ddm_model.R",
  "igt_mod/experimental/igt_mod_pvldecay_tic_model.R",
  "igt_mod/experimental/igt_mod_pvldelta_ddm_b1b_model.R",
  "igt_mod/experimental/igt_mod_pvldelta_ddm_b1p2_model.R",
  "igt_mod/experimental/igt_mod_pvldelta_tic_ddm_model.R",
  "igt_mod/experimental/igt_mod_pvldelta_tic_model.R",
  # working models
  "igt_mod/working/igt_mod_dualp_bbup_ddm_b1p2_model.R",
  "igt_mod/working/igt_mod_pvlboth_ddm_b1b_model.R",
  "igt_mod/working/igt_mod_pvlboth_ddm_b1p2_model.R",
  "igt_mod/working/igt_mod_pvlboth_tic_ddm_model.R",
  "igt_mod/working/igt_mod_pvlboth_tic_model.R"
)

# Function to update a single file
update_file <- function(filepath) {
  cat("Updating:", filepath, "\n")
  
  # Read the file
  content <- readLines(filepath, warn = FALSE)
  full_text <- paste(content, collapse = "\n")
  
  # 1. Update simulate_choices signature if needed
  full_text <- str_replace(
    full_text,
    "simulate_choices = function\\(trials, parameters\\)",
    "simulate_choices = function(trials, parameters, task_params)"
  )
  
  # 2. Update calculate_loglik signature
  # Match the old signature with flexible spacing
  old_sig_pattern <- "calculate_loglik = function\\(trials, choices, RTs, outcomes, parameters(, *\n? *task_params)?\\)"
  new_sig <- "calculate_loglik = function(data, parameters, task_params)"
  
  full_text <- str_replace(
    full_text,
    old_sig_pattern,
    new_sig
  )
  
  # 3. Add variable extraction block after the new signature
  # Find the position right after "calculate_loglik = function(data, parameters, task_params) {"
  extraction_block <- '               # Extract data
               n_trials <- nrow(data)
               choices <- data$choice
               RTs <- data$RT
               outcomes <- data$outcome
               deck_shown <- data$deck_shown
               
               trial_loglik <- numeric(n_trials)'
  
  # Replace the old initialization
  full_text <- str_replace(
    full_text,
    "calculate_loglik = function\\(data, parameters, task_params\\) \\{\n( *)n_trials <- length\\(choices\\)\n( *)trial_loglik <- numeric\\(n_trials\\)",
    paste0("calculate_loglik = function(data, parameters, task_params) {\n", extraction_block)
  )
  
  # 4. Replace trials$deck_shown[t] with deck_shown[t]
  full_text <- str_replace_all(
    full_text,
    "trials\\$deck_shown\\[t\\]",
    "deck_shown[t]"
  )
  
  # Write back to file
  writeLines(full_text, filepath)
  cat("  ✓ Updated successfully\n")
  
  return(TRUE)
}

# Get the current working directory (should be the models directory)
base_dir <- getwd()
cat("Running from:", base_dir, "\n\n")

# Process each file
success_count <- 0
error_count <- 0

for (file in files_to_update) {
  filepath <- file.path(base_dir, file)
  
  if (!file.exists(filepath)) {
    cat("ERROR: File not found:", filepath, "\n")
    error_count <- error_count + 1
    next
  }
  
  tryCatch({
    update_file(filepath)
    success_count <- success_count + 1
  }, error = function(e) {
    cat("ERROR updating", filepath, ":", e$message, "\n")
    error_count <- error_count + 1
  })
}

cat("\n=================================\n")
cat("Update complete!\n")
cat("Successfully updated:", success_count, "files\n")
cat("Errors:", error_count, "files\n")
cat("=================================\n")
