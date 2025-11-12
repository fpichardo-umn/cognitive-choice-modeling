#!/usr/bin/env Rscript

# Model comparison script for mIGT models
# This script compares different models based on their parameter recovery and fit metrics

suppressPackageStartupMessages({
  library(optparse)
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Parse command line arguments
option_list = list(
  make_option(c("-k", "--task"), type="character", default="igt", help="Task name"),
  make_option(c("-m", "--models"), type="character", help="Comma-separated list of models to compare"),
  make_option(c("-g", "--group"), type="character", default="batch_001", help="Group type (sing, hier)"),
  make_option(c("-c", "--cohort"), type="character", default=NULL, help="Cohort identifier"),
  make_option(c("-s", "--session"), type="character", default=NULL, help="Session identifier"),
  make_option(c("-d", "--data_dir"), type="character", default=NULL, help="Directory containing recovery data"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, help="Output directory for comparison"),
  make_option(c("-r", "--render"), action="store_true", default=TRUE, help="Render output to HTML")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("Options used:\n")
dput(opt)

# Source helper functions
source(file.path(here::here(), "scripts", "helpers", "helper_functions_cmdSR.R"))
source(file.path(here::here(), "scripts", "parameter_recovery", "helper_functions_PR.R"))
source(file.path(here::here(), "scripts", "parameter_recovery", "recovery", "recovery.R"))

# Set up directories
dirs <- setup_directories(opt$task)

# Set data directory
data_dir <- if(!is.null(opt$data_dir)) {
  opt$data_dir
} else {
  dirs$REC_SIM_DIR
}

# Set output directory
output_dir <- if(!is.null(opt$output_dir)) {
  opt$output_dir
} else {
  file.path(dirs$REC_SIM_DIR, "comparisons")
}

# Create output directory if it doesn't exist
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Parse models
models <- unlist(strsplit(opt$models, ","))
cat("Comparing models:", paste(models, collapse=", "), "\n")

# Function to load recovery data
load_recovery_data <- function(model, task, group, cohort, session, data_dir) {
  csv_file <- file.path(
    data_dir,
    generate_bids_filename(
      prefix = NULL,
      task = task,
      group = group,
      model = model,
      cohort = cohort,
      ses = session,
      additional_tags = list(
        "type" = "rec",
        "desc" = "data"
      ),
      ext = "csv"
    )
  )
  
  if(!file.exists(csv_file)) {
    warning("Recovery data file not found: ", csv_file)
    return(NULL)
  }
  
  data <- read.csv(csv_file)
  data$model <- model  # Add model identifier
  return(data)
}

# Load recovery data for all models
recovery_data_list <- list()
for(model in models) {
  data <- load_recovery_data(model, opt$task, opt$group, opt$cohort, opt$session, data_dir)
  if(!is.null(data)) {
    recovery_data_list[[model]] <- data
  }
}

if(length(recovery_data_list) == 0) {
  stop("No recovery data found for any models")
}

# Combine all recovery data
all_recovery_data <- bind_rows(recovery_data_list)

# Calculate recovery statistics for each model
recovery_stats <- list()
for(model in models) {
  model_data <- all_recovery_data[all_recovery_data$model == model, ]
  if(nrow(model_data) > 0) {
    recovery_stats[[model]] <- calculate_recovery_statistics(model_data)
  }
}

# Create comparison RMD file
rmd_file <- file.path(
  output_dir,
  generate_bids_filename(
    prefix = NULL,
    task = opt$task,
    group = opt$group,
    model = "comparison",
    cohort = opt$cohort,
    ses = opt$session,
    additional_tags = list(
      "type" = "rec",
      "desc" = "comparison"
    ),
    ext = "Rmd"
  )
)

# Get template path
template_dir <- file.path(dirs$PR_DIR, "recovery", "templates")
template_path <- file.path(template_dir, "model_comparison_template.Rmd")

# If template doesn't exist, need to create it first using the existing scripts as a guide
if(!file.exists(template_path)) {
  warning("Model comparison template not found, creating it...")
  # Create model comparison template - simplified example content
  template_content <- '---
title: "Model Comparison: {{TASK}} - {{GROUP}}"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: cosmo
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
library(reshape2)
library(here)

# Data for comparison
models <- c({{MODELS}})
task <- "{{TASK}}"
group <- "{{GROUP}}"
recovery_data <- read.csv("{{DATA_FILE}}")
```

## Model Comparison for {{TASK}} - {{GROUP}}

This analysis compares different computational models applied to the {{TASK}} task with group type {{GROUP}}.

### Overall Recovery Performance

```{r overall_comparison, fig.width=10, fig.height=8}
# Summarize standardized metrics by model
model_summary <- recovery_data %>%
  group_by(model) %>%
  summarize(
    correlation = cor(true_value, recovered_value, use="pairwise.complete.obs"),
    rmse = sqrt(mean((true_value - recovered_value)^2, na.rm=TRUE)),
    bias = mean(recovered_value - true_value, na.rm=TRUE),
    relative_bias = mean((recovered_value - true_value) / pmax(0.0001, abs(true_value)), na.rm=TRUE)
  ) %>%
  arrange(desc(correlation))

# Display overall metrics
kable(model_summary, caption = "Overall Model Recovery Metrics") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))

# Plot standardized metrics comparison
model_summary_long <- model_summary %>%
  pivot_longer(cols = -model, names_to = "metric", values_to = "value")

ggplot(model_summary_long, aes(x = model, y = value, fill = model)) +
  geom_bar(stat = "identity") +
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Model Comparison - Recovery Metrics", x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

### Parameter-Specific Recovery

```{r parameter_comparison, fig.width=12, fig.height=10}
# Compare parameter recovery by parameter type
param_recovery <- recovery_data %>%
  group_by(model, parameter) %>%
  summarize(
    correlation = cor(true_value, recovered_value, use="pairwise.complete.obs"),
    rmse = sqrt(mean((true_value - recovered_value)^2, na.rm=TRUE)),
    bias = mean(recovered_value - true_value, na.rm=TRUE),
    .groups = "drop"
  )

# Display parameter metrics
kable(param_recovery, caption = "Parameter Recovery by Model") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  collapse_rows(columns = 1, valign = "top")

# Plot parameter recovery metrics
param_recovery_long <- param_recovery %>%
  pivot_longer(cols = c(correlation, rmse, bias), 
               names_to = "metric", values_to = "value")

ggplot(param_recovery_long, aes(x = parameter, y = value, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "Parameter Recovery by Model", x = "Parameter", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

### True vs. Recovered Parameter Plots

```{r parameter_plots, fig.width=12, fig.height=10}
# Plot true vs. recovered values by model and parameter
# Focus on common parameters across models
common_params <- recovery_data %>%
  group_by(parameter) %>%
  summarize(num_models = length(unique(model))) %>%
  filter(num_models > 1) %>%
  pull(parameter)

if(length(common_params) > 0) {
  recovery_data %>%
    filter(parameter %in% common_params) %>%
    ggplot(aes(x = true_value, y = recovered_value, color = model)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    facet_wrap(~parameter, scales = "free") +
    labs(title = "Parameter Recovery Comparison", 
         x = "True Parameter Value", 
         y = "Recovered Parameter Value") +
    theme_minimal()
}
```

### Model Decision Characteristics

```{r decision_characteristics, fig.width=10, fig.height=8}
# Here you would include model-specific behavioral predictions
# This would require additional data about model simulations
```

### Conclusion

This analysis compares the parameter recovery performance of different computational models for the {{TASK}} task. Key metrics include correlation between true and recovered parameters, RMSE, and bias.

The comparison helps identify which models have better parameter recovery, which is an important aspect of model evaluation along with goodness-of-fit and theoretical plausibility.
'

  # Write template file
  writeLines(template_content, template_path)
}

# Read the template
template_text <- readLines(template_path, warn = FALSE)
template_text <- paste(template_text, collapse = "\n")

# Create combined data file for all models
combined_data_file <- file.path(
  output_dir,
  generate_bids_filename(
    prefix = NULL,
    task = opt$task,
    group = opt$group,
    model = "all",
    cohort = opt$cohort,
    ses = opt$session,
    additional_tags = list(
      "type" = "rec",
      "desc" = "data"
    ),
    ext = "csv"
  )
)

# Write combined data to CSV
write.csv(all_recovery_data, combined_data_file, row.names = FALSE)

# Replace placeholders in template
rmd_content <- gsub("\\{\\{TASK\\}\\}", opt$task, template_text)
rmd_content <- gsub("\\{\\{GROUP\\}\\}", opt$group, rmd_content)
rmd_content <- gsub("\\{\\{MODELS\\}\\}", paste(shQuote(models), collapse = ", "), rmd_content)
rmd_content <- gsub("\\{\\{DATA_FILE\\}\\}", basename(combined_data_file), rmd_content)

# Write to output file
writeLines(rmd_content, rmd_file)

# Render to HTML if requested
if(opt$render) {
  html_file <- gsub("\\.Rmd$", ".html", rmd_file)
  cat("Rendering HTML report:", html_file, "\n")
  rmarkdown::render(rmd_file, output_file = html_file)
}

cat("\nModel comparison complete.\n")
cat("- Comparison RMD file:", rmd_file, "\n")
if(opt$render) {
  cat("- Comparison HTML file:", html_file, "\n")
}
