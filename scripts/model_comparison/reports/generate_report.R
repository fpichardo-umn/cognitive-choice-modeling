#!/usr/bin/env Rscript

#' Generate Model Comparison Report
#' @description Create comprehensive HTML report for model comparison analysis

suppressPackageStartupMessages({
  library(here)
  library(rmarkdown)
  library(knitr)
  library(dplyr)
})

# Load helper functions
source(file.path(here::here(), "scripts", "model_comparison", "helpers", "model_comparison_helpers.R"))

#' Generate comprehensive model comparison report
#' @param analysis_results Results from all analysis modules
#' @param comparison_data Original comparison data
#' @param models_by_type Models organized by type
#' @param output_dir Directory for report output
#' @param results_file Path to saved results RDS file
#' @param task Task name
#' @param cohort Cohort name
#' @param session Session name (optional)
#' @param comparison_name Name of this comparison
#' @param render_html Whether to render HTML
#' @return List with report file paths
generate_model_comparison_report <- function(analysis_results, comparison_data, models_by_type, 
                                           output_dir, results_file, task, cohort, session = NULL, 
                                           comparison_name, render_html = TRUE) {
  message("Generating model comparison report...")
  
  # Create report filename
  rmd_file_base <- generate_bids_filename(
    prefix = NULL,
    task = task,
    cohort = cohort,
    ses = session,
    group = "comparison",
    model = NULL,
    additional_tags = list(
      comparison = comparison_name,
      type = "modelcomp",
      desc = "integrated"
    ),
    ext = "Rmd"
  )
  
  rmd_file <- file.path(output_dir, rmd_file_base)
  html_file <- gsub("Rmd", "html", rmd_file_base)
  
  # Generate RMD content
  rmd_content <- create_report_rmd_content(
    analysis_results, comparison_data, models_by_type,
    task, cohort, session, comparison_name, results_file
  )
  
  # Write RMD file
  writeLines(rmd_content, rmd_file)
  
  # Render to HTML if requested
  if (render_html) {
    message("Rendering HTML report...")
    
    tryCatch({
      rmarkdown::render(
        input = rmd_file,
        output_file = html_file,
        quiet = TRUE,
        envir = new.env()
      )
      message("HTML report generated: ", html_file)
    }, error = function(e) {
      warning("Failed to render HTML report: ", e$message)
      html_file <- NULL
    })
  }
  
  return(list(
    rmd_file = rmd_file,
    html_file = html_file
  ))
}

#' Create RMD content for the report
#' @param analysis_results Results from all analysis modules
#' @param comparison_data Original comparison data
#' @param models_by_type Models organized by type
#' @param task Task name
#' @param cohort Cohort name
#' @param session Session name
#' @param comparison_name Comparison name
#' @param results_file Path to saved results RDS file
#' @return Character vector with RMD content
create_report_rmd_content <- function(analysis_results, comparison_data, models_by_type,
                                    task, cohort, session, comparison_name, results_file) {
  
  # Load template components
  template_dir <- file.path(here::here(), "scripts", "model_comparison", "reports", "templates")
  
  # Use main template and load section templates
  main_template <- read_template_section(template_dir, "main_comparison_template.Rmd")
  ic_section <- read_template_section(template_dir, "ic_section_template.Rmd")
  recovery_section <- read_template_section(template_dir, "recovery_section_template.Rmd")
  ppc_section <- read_template_section(template_dir, "ppc_section_template.Rmd")
  model_type_section <- read_template_section(template_dir, "model_type_section_template.Rmd")
  task_specific_section <- read_template_section(template_dir, "task_specific_section_template.Rmd")
  
  # Use the absolute path to results file (passed as parameter)
  results_file_path <- results_file
  
  # Replace placeholders in templates
  placeholders <- list(
    TASK = task,
    COHORT = cohort,
    SESSION = session %||% "none",
    COMPARISON_NAME = comparison_name,
    TIMESTAMP = as.character(Sys.time()),
    N_MODELS = length(comparison_data),
    RESULTS_FILE = results_file_path,
    IC_SECTION = paste(ic_section, collapse = "\n"),
    RECOVERY_SECTION = paste(recovery_section, collapse = "\n"),
    PPC_SECTION = paste(ppc_section, collapse = "\n"),
    MODEL_TYPE_SECTION = paste(model_type_section, collapse = "\n"),
    TASK_SPECIFIC_SECTION = paste(task_specific_section, collapse = "\n"),
    MODEL_PROFILES_SECTION = create_model_profiles_section(comparison_data)
  )
  
  # Process main template with all placeholders
  full_content <- replace_placeholders(main_template, placeholders)
  
  return(paste(full_content, collapse = "\n"))
}

#' Read template section from file
#' @param template_dir Template directory
#' @param filename Template filename
#' @return Character vector with template content
read_template_section <- function(template_dir, filename) {
  template_path <- file.path(template_dir, filename)
  
  if (file.exists(template_path)) {
    return(readLines(template_path, warn = FALSE))
  } else {
    # Return default content if template not found
    warning("Template not found: ", template_path, ". Using default content.")
    return(get_default_template_content(filename))
  }
}

#' Replace placeholders in template content
#' @param template_lines Character vector of template lines
#' @param placeholders Named list of placeholder values
#' @return Character vector with placeholders replaced
replace_placeholders <- function(template_lines, placeholders) {
  content <- paste(template_lines, collapse = "\n")
  
  for (placeholder_name in names(placeholders)) {
    placeholder_pattern <- paste0("\\{\\{", placeholder_name, "\\}\\}")
    content <- gsub(placeholder_pattern, placeholders[[placeholder_name]], content)
  }
  
  return(content)
}

#' Create model profiles section dynamically
#' @param comparison_data Original comparison data
#' @return Character vector with model profiles section
create_model_profiles_section <- function(comparison_data) {
  profiles_content <- c(
    "```{r model-profiles}",
    "# Create detailed profile for each model",
    "model_names <- names(comparison_data)",
    "",
    "for (model_name in model_names) {",
    "  cat('\n## Model:', model_name, '\n\n')",
    "  ",
    "  # Model metadata",
    "  model_type <- classify_model_type(model_name)",
    "  cat('**Type:** ', model_type, '\n\n')",
    "  ",
    "  # IC performance",
    "  if ('ic' %in% names(analysis_results) && nrow(analysis_results$ic$overall_ranking) > 0) {",
    "    ic_info <- analysis_results$ic$overall_ranking %>% filter(model == model_name)",
    "    if (nrow(ic_info) > 0) {",
    "      cat('**Information Criteria:**\n')",
    "      cat('- Rank:', ic_info$rank, '\n')",
    "      cat('- Performance Tier:', ic_info$performance_tier, '\n')",
    "      cat('- Model Weight:', round(ic_info$model_weight, 3), '\n\n')",
    "    }",
    "  }",
    "  ",
    "  # Recovery performance", 
    "  if ('recovery' %in% names(analysis_results) && nrow(analysis_results$recovery$model_summary) > 0) {",
    "    recovery_info <- analysis_results$recovery$model_summary %>% filter(model == model_name)",
    "    if (nrow(recovery_info) > 0) {",
    "      cat('**Parameter Recovery:**\n')",
    "      cat('- Mean Correlation:', round(recovery_info$mean_correlation, 3), '\n')",
    "      cat('- Recovery Quality:', recovery_info$recovery_quality, '\n\n')",
    "    }",
    "  }",
    "  ",
    "  # PPC performance",
    "  if ('ppc' %in% names(analysis_results) && nrow(analysis_results$ppc$model_summary) > 0) {",
    "    ppc_info <- analysis_results$ppc$model_summary %>% filter(model == model_name)",
    "    if (nrow(ppc_info) > 0) {",
    "      cat('**PPC Performance:**\n')",
    "      cat('- Mean PPP:', round(ppc_info$overall_mean_ppp, 3), '\n')",
    "      cat('- Proportion Extreme:', round(ppc_info$overall_proportion_extreme, 3), '\n')",
    "      cat('- PPC Quality:', ppc_info$model_quality, '\n')",
    "      cat('- Best Domain:', ppc_info$best_domain, '\n')",
    "      cat('- Worst Domain:', ppc_info$worst_domain, '\n\n')",
    "    }",
    "  }",
    "  ",
    "  cat('---\n')",
    "}",
    "```"
  )
  
  return(paste(profiles_content, collapse = "\n"))
}

#' Get default template content for missing templates
#' @param filename Template filename
#' @return Default template content
get_default_template_content <- function(filename) {
  switch(filename,
    "main_comparison_template.Rmd" = get_default_main_template(),
    "ic_section_template.Rmd" = get_default_ic_analysis_template(),
    "recovery_section_template.Rmd" = get_default_recovery_analysis_template(),
    "ppc_section_template.Rmd" = get_default_ppc_analysis_template(),
    "model_type_section_template.Rmd" = get_default_model_type_template(),
    "task_specific_section_template.Rmd" = get_default_task_specific_template(),
    paste("# Unknown Template:", filename)
  )
}

#' Default main template
get_default_main_template <- function() {
  c(
    "---",
    "title: 'Model Comparison Report: {{TASK}} - {{COHORT}}'",
    "subtitle: '{{COMPARISON_NAME}}'",
    "date: '`r Sys.Date()`'",
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    toc_depth: 3",
    "    theme: flatly",
    "    code_folding: hide",
    "    fig_width: 10",
    "    fig_height: 6",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(dplyr)",
    "library(ggplot2)",
    "library(knitr)",
    "library(DT)",
    "library(here)",
    "",
    "# Load analysis results",
    "results_file <- '{{RESULTS_FILE}}'",
    "all_results <- readRDS(results_file)",
    "comparison_data <- all_results$comparison_data",
    "analysis_results <- all_results$analysis_results",
    "models_by_type <- all_results$models_by_type",
    "config <- all_results$config",
    "",
    "# Load helper functions",
    "source(file.path(here::here(), 'scripts', 'model_comparison', 'helpers', 'model_comparison_helpers.R'))",
    "```",
    "",
    "# Executive Summary",
    "",
    "{{IC_SECTION}}",
    "",
    "{{RECOVERY_SECTION}}",
    "",
    "{{PPC_SECTION}}",
    "",
    "{{MODEL_TYPE_SECTION}}",
    "",
    "{{TASK_SPECIFIC_SECTION}}",
    "",
    "# Detailed Model Profiles",
    "",
    "{{MODEL_PROFILES_SECTION}}"
  )
}

#' Default model type template
get_default_model_type_template <- function() {
  c(
    "# Model Type Comparisons",
    "",
    "Analysis of performance differences across model types (RL, SSM, hybrid).",
    "",
    "```{r model-type-analysis}",
    "if (length(models_by_type) > 1) {",
    "  cat('Multiple model types available for comparison.\n')",
    "} else {",
    "  cat('Only one model type in this comparison.\n')",
    "}",
    "```"
  )
}

#' Default task specific template
get_default_task_specific_template <- function() {
  c(
    "# Task-Specific Analysis",
    "",
    "Analysis tailored to the specific task requirements.",
    "",
    "```{r task-specific-analysis}",
    "cat('Task:', config$task, '\n')",
    "```"
  )
}

#' Default header template
get_default_header_template <- function() {
  c(
    "---",
    "title: 'Model Comparison Report: {{TASK}} - {{COHORT}}'",
    "author: 'Model Comparison Pipeline'",
    "date: '{{TIMESTAMP}}'",
    "output:",
    "  html_document:",
    "    theme: flatly",
    "    toc: true",
    "    toc_float: true",
    "    code_folding: hide",
    "    fig_width: 10",
    "    fig_height: 6",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(dplyr)",
    "library(ggplot2)",
    "library(knitr)",
    "library(DT)",
    "",
    "# Load analysis results",
    "results_file <- file.path(dirname(getwd()), 'data', 'model_comparison_results.rds')",
    "if (file.exists(results_file)) {",
    "  all_results <- readRDS(results_file)",
    "  analysis_results <- all_results$analysis_results",
    "  comparison_data <- all_results$comparison_data",
    "  models_by_type <- all_results$models_by_type",
    "} else {",
    "  stop('Analysis results file not found')",
    "}",
    "```"
  )
}

#' Default executive summary template
get_default_executive_summary_template <- function() {
  c(
    "# Executive Summary",
    "",
    "This report presents a comprehensive comparison of {{N_MODELS}} computational models for the **{{TASK}}** task using data from the **{{COHORT}}** cohort.",
    "",
    "```{r executive-summary}",
    "# Create executive summary table",
    "if ('ic' %in% names(analysis_results) && nrow(analysis_results$ic$overall_ranking) > 0) {",
    "  # Top 5 models by information criteria",
    "  top_models <- head(analysis_results$ic$overall_ranking, 5)",
    "  ",
    "  # Get recovery and PPC info for top models",
    "  if ('recovery' %in% names(analysis_results)) {",
    "    recovery_info <- analysis_results$recovery$model_summary %>%",
    "      select(model, recovery_quality = recovery_quality, mean_correlation)",
    "    top_models <- top_models %>% left_join(recovery_info, by = 'model')",
    "  }",
    "  ",
    "  if ('ppc' %in% names(analysis_results)) {",
    "    ppc_info <- analysis_results$ppc$model_summary %>%",
    "      select(model, ppc_quality = model_quality, overall_proportion_extreme)",
    "    top_models <- top_models %>% left_join(ppc_info, by = 'model')",
    "  }",
    "  ",
    "  # Display summary table",
    "  kable(top_models %>% select(Rank = rank, Model = model, Type = model_type, ",
    "                             everything(), -performance_tier),",
    "        caption = 'Top 5 Models Summary', digits = 3)",
    "} else {",
    "  cat('No information criteria results available.')",
    "}",
    "```",
    "",
    "## Key Findings",
    "",
    "```{r key-findings}",
    "# Generate key findings",
    "findings <- list()",
    "",
    "# Best model overall",
    "if ('ic' %in% names(analysis_results) && nrow(analysis_results$ic$overall_ranking) > 0) {",
    "  best_model <- analysis_results$ic$overall_ranking$model[1]",
    "  findings$best_model <- paste('**Best Model (IC):**', best_model)",
    "}",
    "",
    "# Parameter recovery insights",
    "if ('recovery' %in% names(analysis_results) && nrow(analysis_results$recovery$group_summary) > 0) {",
    "  best_group <- analysis_results$recovery$group_summary$group[1]",
    "  worst_group <- tail(analysis_results$recovery$group_summary$group, 1)",
    "  findings$recovery <- paste('**Parameter Recovery:** Best recovered group:', best_group, '| Worst:', worst_group)",
    "}",
    "",
    "# PPC insights",
    "if ('ppc' %in% names(analysis_results) && nrow(analysis_results$ppc$domain_summary) > 0) {",
    "  hardest_domain <- analysis_results$ppc$domain_summary %>%",
    "    arrange(desc(proportion_extreme)) %>%",
    "    slice(1) %>%",
    "    pull(domain)",
    "  findings$ppc <- paste('**PPC Performance:** Most challenging behavioral domain:', hardest_domain)",
    "}",
    "",
    "# Display findings",
    "for (finding in findings) {",
    "  cat(finding, '\\n\\n')",
    "}",
    "```"
  )
}

#' Default IC analysis template
get_default_ic_analysis_template <- function() {
  c(
    "# Information Criteria Analysis",
    "",
    "Model selection based on information criteria (LOOIC/WAIC) provides an objective measure of model fit while penalizing complexity.",
    "",
    "## Overall Model Ranking",
    "",
    "```{r ic-ranking-table}",
    "if ('ic' %in% names(analysis_results) && nrow(analysis_results$ic$overall_ranking) > 0) {",
    "  ranking_table <- analysis_results$ic$overall_ranking %>%",
    "    select(Rank = rank, Model = model, Type = model_type, ",
    "           matches('estimate$'), matches('^delta_'), Weight = model_weight, Tier = performance_tier)",
    "  ",
    "  DT::datatable(ranking_table, options = list(pageLength = 15)) %>%",
    "    DT::formatRound(columns = c(4:6), digits = 2)",
    "} else {",
    "  cat('No information criteria results available.')",
    "}",
    "```",
    "",
    "## Model Performance Tiers",
    "",
    "```{r ic-tiers}",
    "if ('ic' %in% names(analysis_results) && nrow(analysis_results$ic$overall_ranking) > 0) {",
    "  tier_summary <- analysis_results$ic$overall_ranking %>%",
    "    count(performance_tier, name = 'n_models') %>%",
    "    arrange(match(performance_tier, c('top_tier', 'competitive', 'considerable_difference', 'strong_difference', 'decisive_difference')))",
    "  ",
    "  kable(tier_summary, col.names = c('Performance Tier', 'Number of Models'),",
    "        caption = 'Model Performance Distribution')",
    "}",
    "```",
    "",
    "## Rankings by Model Type",
    "",
    "```{r ic-by-type}",
    "if ('ic' %in% names(analysis_results) && length(analysis_results$ic$by_type_ranking) > 0) {",
    "  for (type_name in names(analysis_results$ic$by_type_ranking)) {",
    "    cat('\\n### ', toupper(type_name), ' Models\\n\\n')",
    "    ",
    "    type_ranking <- analysis_results$ic$by_type_ranking[[type_name]] %>%",
    "      select(Rank = rank, Model = model, matches('estimate$'), matches('^delta_'), Weight = model_weight)",
    "    ",
    "    print(kable(type_ranking, caption = paste(toupper(type_name), 'Model Rankings'), digits = 3))",
    "  }",
    "} else {",
    "  cat('No model type rankings available.')",
    "}",
    "```"
  )
}

#' Default recovery analysis template
get_default_recovery_analysis_template <- function() {
  c(
    "# Parameter Recovery Analysis",
    "",
    "Parameter recovery analysis tests whether model parameters can be accurately recovered from simulated data, organized by psychological construct groups.",
    "",
    "## Recovery Quality by Parameter Group",
    "",
    "```{r recovery-group-summary}",
    "if ('recovery' %in% names(analysis_results) && nrow(analysis_results$recovery$group_summary) > 0) {",
    "  group_table <- analysis_results$recovery$group_summary %>%",
    "    select(Group = group, 'N Models' = n_models, 'Mean Correlation' = mean_correlation,",
    "           'Min Correlation' = min_correlation, 'Max Correlation' = max_correlation,",
    "           Quality = recovery_quality)",
    "  ",
    "  DT::datatable(group_table, options = list(pageLength = 10)) %>%",
    "    DT::formatRound(columns = 3:5, digits = 3)",
    "} else {",
    "  cat('No parameter recovery results available.')",
    "}",
    "```",
    "",
    "## Recovery Quality by Model",
    "",
    "```{r recovery-model-summary}",
    "if ('recovery' %in% names(analysis_results) && nrow(analysis_results$recovery$model_summary) > 0) {",
    "  model_table <- analysis_results$recovery$model_summary %>%",
    "    select(Model = model, Type = model_type, 'N Parameters' = n_parameters,", 
    "           'Mean Correlation' = mean_correlation, 'Min Correlation' = min_correlation,",
    "           'Max Correlation' = max_correlation, Quality = recovery_quality)",
    "  ",
    "  DT::datatable(model_table, options = list(pageLength = 15)) %>%",
    "    DT::formatRound(columns = 4:6, digits = 3)",
    "} else {",
    "  cat('No model-level recovery results available.')",
    "}",
    "```",
    "",
    "## Cross-Model Type Comparison",
    "",
    "```{r recovery-cross-type}",
    "if ('recovery' %in% names(analysis_results) && nrow(analysis_results$recovery$cross_model_comparison) > 0) {",
    "  cross_type_table <- analysis_results$recovery$cross_model_comparison %>%",
    "    select('Model Type' = model_type, Group = group, 'N Models' = n_models,",
    "           'Mean Correlation' = mean_correlation, 'SD Correlation' = sd_correlation,",
    "           Quality = recovery_quality)",
    "  ",
    "  DT::datatable(cross_type_table, options = list(pageLength = 15)) %>%",
    "    DT::formatRound(columns = 4:5, digits = 3)",
    "} else {",
    "  cat('No cross-model type recovery comparison available.')",
    "}",
    "```",
    "",
    "## Key Recovery Insights",
    "",
    "```{r recovery-insights}",
    "if ('recovery' %in% names(analysis_results)) {",
    "  insights <- list()",
    "  ",
    "  # Best and worst groups",
    "  if (length(analysis_results$recovery$best_worst_groups) > 0) {",
    "    if (nrow(analysis_results$recovery$best_worst_groups$best_groups) > 0) {",
    "      best_groups <- head(analysis_results$recovery$best_worst_groups$best_groups$group, 3)",
    "      insights$best <- paste('**Best Recovered Groups:**', paste(best_groups, collapse = ', '))",
    "    }",
    "    ",
    "    if (nrow(analysis_results$recovery$best_worst_groups$worst_groups) > 0) {",
    "      worst_groups <- head(analysis_results$recovery$best_worst_groups$worst_groups$group, 3)",
    "      insights$worst <- paste('**Worst Recovered Groups:**', paste(worst_groups, collapse = ', '))",
    "    }",
    "  }",
    "  ",
    "  # Display insights",
    "  for (insight in insights) {",
    "    cat(insight, '\\n\\n')",
    "  }",
    "}",
    "```"
  )
}

#' Default PPC analysis template
get_default_ppc_analysis_template <- function() {
  c(
    "# Posterior Predictive Checks Analysis",
    "",
    "Posterior predictive checks assess how well models reproduce observed behavioral patterns, organized by behavioral domains.",
    "",
    "## PPC Performance by Behavioral Domain",
    "",
    "```{r ppc-domain-summary}",
    "if ('ppc' %in% names(analysis_results) && nrow(analysis_results$ppc$domain_summary) > 0) {",
    "  domain_table <- analysis_results$ppc$domain_summary %>%",
    "    select(Domain = domain, 'N Models' = n_models, 'Mean PPP' = mean_ppp,",
    "           'Proportion Extreme' = proportion_extreme, 'N Statistics' = n_statistics,",
    "           Quality = domain_quality)",
    "  ",
    "  DT::datatable(domain_table, options = list(pageLength = 10)) %>%",
    "    DT::formatRound(columns = 3:4, digits = 3)",
    "} else {",
    "  cat('No PPC domain results available.')",
    "}",
    "```",
    "",
    "## PPC Performance by Model",
    "",
    "```{r ppc-model-summary}",
    "if ('ppc' %in% names(analysis_results) && nrow(analysis_results$ppc$model_summary) > 0) {",
    "  model_table <- analysis_results$ppc$model_summary %>%",
    "    select(Model = model, Type = model_type, 'N Domains' = n_domains,",
    "           'Mean PPP' = overall_mean_ppp, 'Proportion Extreme' = overall_proportion_extreme,",
    "           Quality = model_quality, 'Worst Domain' = worst_domain, 'Best Domain' = best_domain)",
    "  ",
    "  DT::datatable(model_table, options = list(pageLength = 15)) %>%",
    "    DT::formatRound(columns = 4:5, digits = 3)",
    "} else {",
    "  cat('No model-level PPC results available.')",
    "}",
    "```",
    "",
    "## Extreme PPC Failures",
    "",
    "```{r ppc-extreme-failures}",
    "if ('ppc' %in% names(analysis_results) && nrow(analysis_results$ppc$extreme_failures) > 0) {",
    "  failure_summary <- analysis_results$ppc$extreme_failures %>%",
    "    count(model, model_type, failure_severity, name = 'n_failures') %>%",
    "    pivot_wider(names_from = failure_severity, values_from = n_failures, values_fill = 0) %>%",
    "    arrange(desc(severe), desc(moderate))",
    "  ",
    "  DT::datatable(failure_summary, caption = 'Extreme PPC Failures by Model',",
    "                options = list(pageLength = 10))",
    "} else {",
    "  cat('No extreme PPC failures detected.')",
    "}",
    "```",
    "",
    "## Behavioral Pattern Analysis",
    "",
    "```{r ppc-behavioral-patterns}",
    "if ('ppc' %in% names(analysis_results) && 'behavioral_patterns' %in% names(analysis_results$ppc)) {",
    "  patterns <- analysis_results$ppc$behavioral_patterns",
    "  ",
    "  if ('difficult_domains' %in% names(patterns) && nrow(patterns$difficult_domains) > 0) {",
    "    cat('### Most Challenging Behavioral Domains\\n\\n')",
    "    ",
    "    difficult_table <- patterns$difficult_domains %>%",
    "      select('Model Type' = model_type, Domain = domain, 'N Models' = n_models,",
    "             'Mean Extreme Proportion' = mean_proportion_extreme) %>%",
    "      arrange(desc(`Mean Extreme Proportion`))",
    "    ",
    "    print(kable(difficult_table, digits = 3, caption = 'Challenging Domains by Model Type'))",
    "  }",
    "}",
    "```"
  )
}

#' Default detailed profiles template
get_default_detailed_profiles_template <- function() {
  c(
    "# Detailed Model Profiles",
    "",
    "Individual profiles for each model showing comprehensive performance across all analysis dimensions.",
    "",
    "```{r model-profiles}",
    "# Create detailed profile for each model",
    "model_names <- names(comparison_data)",
    "",
    "for (model_name in model_names) {",
    "  cat('\\n## Model:', model_name, '\\n\\n')",
    "  ",
    "  # Model metadata",
    "  model_type <- classify_model_type(model_name)",
    "  cat('**Type:** ', model_type, '\\n\\n')",
    "  ",
    "  # IC performance",
    "  if ('ic' %in% names(analysis_results) && nrow(analysis_results$ic$overall_ranking) > 0) {",
    "    ic_info <- analysis_results$ic$overall_ranking %>% filter(model == model_name)",
    "    if (nrow(ic_info) > 0) {",
    "      cat('**Information Criteria:**\\n')",
    "      cat('- Rank:', ic_info$rank, '\\n')",
    "      cat('- Performance Tier:', ic_info$performance_tier, '\\n')",
    "      cat('- Model Weight:', round(ic_info$model_weight, 3), '\\n\\n')",
    "    }",
    "  }",
    "  ",
    "  # Recovery performance", 
    "  if ('recovery' %in% names(analysis_results) && nrow(analysis_results$recovery$model_summary) > 0) {",
    "    recovery_info <- analysis_results$recovery$model_summary %>% filter(model == model_name)",
    "    if (nrow(recovery_info) > 0) {",
    "      cat('**Parameter Recovery:**\\n')",
    "      cat('- Mean Correlation:', round(recovery_info$mean_correlation, 3), '\\n')",
    "      cat('- Recovery Quality:', recovery_info$recovery_quality, '\\n\\n')",
    "    }",
    "  }",
    "  ",
    "  # PPC performance",
    "  if ('ppc' %in% names(analysis_results) && nrow(analysis_results$ppc$model_summary) > 0) {",
    "    ppc_info <- analysis_results$ppc$model_summary %>% filter(model == model_name)",
    "    if (nrow(ppc_info) > 0) {",
    "      cat('**PPC Performance:**\\n')",
    "      cat('- Mean PPP:', round(ppc_info$overall_mean_ppp, 3), '\\n')",
    "      cat('- Proportion Extreme:', round(ppc_info$overall_proportion_extreme, 3), '\\n')",
    "      cat('- PPC Quality:', ppc_info$model_quality, '\\n')",
    "      cat('- Best Domain:', ppc_info$best_domain, '\\n')",
    "      cat('- Worst Domain:', ppc_info$worst_domain, '\\n\\n')",
    "    }",
    "  }",
    "  ",
    "  cat('---\\n')",
    "}",
    "```"
  )
}

#' Default technical appendix template
get_default_technical_appendix_template <- function() {
  c(
    "# Technical Appendix",
    "",
    "## Analysis Configuration",
    "",
    "```{r analysis-config}",
    "config_info <- list(",
    "  'Task' = '{{TASK}}',",
    "  'Cohort' = '{{COHORT}}',",
    "  'Session' = '{{SESSION}}',",
    "  'Comparison Name' = '{{COMPARISON_NAME}}',",
    "  'Analysis Timestamp' = '{{TIMESTAMP}}',",
    "  'Number of Models' = {{N_MODELS}}",
    ")",
    "",
    "kable(data.frame(Parameter = names(config_info), Value = unlist(config_info)),",
    "      caption = 'Analysis Configuration')",
    "```",
    "",
    "## Data Sources",
    "",
    "```{r data-sources}",
    "# Information about data sources",
    "cat('### Parameter Recovery Data\\n')",
    "cat('- Source: Outputs/{{TASK}}/validation/parameter_recovery/analysis/\\n')",
    "cat('- Format: CSV files with true vs recovered parameter values\\n\\n')",
    "",
    "cat('### PPC Data\\n')",
    "cat('- Source: Outputs/{{TASK}}/validation/ppc/cohort-{{COHORT}}/ses-{{SESSION}}/stats/\\n')",
    "cat('- Format: CSV files with PPP statistics\\n\\n')",
    "",
    "cat('### Information Criteria Data\\n')",
    "cat('- Source: Outputs/{{TASK}}/validation/ppc/cohort-{{COHORT}}/ses-{{SESSION}}/loglik/\\n')",
    "cat('- Format: RDS files with LOOIC/WAIC estimates\\n\\n')",
    "",
    "cat('### Model Comparison Results\\n')",
    "cat('- Source: Outputs/{{TASK}}/model_comparison/cohort-{{COHORT}}/ses-{{SESSION}}/{{COMPARISON_NAME}}/data/\\n')",
    "cat('- Format: RDS file with consolidated analysis results\\n\\n')",
    "```",
    "",
    "## Parameter Group Definitions",
    "",
    "```{r parameter-groups}",
    "# Load and display parameter groups",
    "config <- load_parameter_groups_config()",
    "groups_df <- data.frame()",
    "",
    "for (group_name in names(config$parameter_groups)) {",
    "  group_info <- config$parameter_groups[[group_name]]",
    "  group_df <- data.frame(",
    "    Group = group_name,",
    "    Description = group_info$description,",
    "    Parameters = paste(group_info$parameters, collapse = ', ')",
    "  )",
    "  groups_df <- rbind(groups_df, group_df)",
    "}",
    "",
    "DT::datatable(groups_df, caption = 'Parameter Group Definitions',",
    "              options = list(pageLength = 10))",
    "```",
    "",
    "## Analysis Methods",
    "",
    "### Parameter Recovery",
    "- **Metric**: Pearson correlation between true and recovered parameters",
    "- **Grouping**: Parameters organized by psychological construct",
    "- **Quality Thresholds**: Excellent (≥0.8), Good (≥0.6), Acceptable (≥0.4), Poor (<0.4)",
    "",
    "### Posterior Predictive Checks",
    "- **Metric**: Posterior Predictive P-values (PPP)",
    "- **Extreme Threshold**: PPP < 0.05 or PPP > 0.95",
    "- **Domains**: Behavioral patterns grouped by psychological relevance",
    "",
    "### Information Criteria",
    "- **Metric**: LOOIC (Leave-One-Out Information Criterion)",
    "- **Model Weights**: Calculated using Akaike weights",
    "- **Performance Tiers**: Based on ΔLOOIC relative to best model",
    "",
    "## Session Information",
    "",
    "```{r session-info}",
    "sessionInfo()",
    "```"
  )
}
