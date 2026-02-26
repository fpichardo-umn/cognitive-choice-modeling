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
    PPC_SECTION = paste(ppc_section, collapse = "\n")
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

#' Get default template content for missing templates
#' @param filename Template filename
#' @return Default template content
get_default_template_content <- function(filename) {
  switch(filename,
    "main_comparison_template.Rmd" = get_default_main_template(),
    "ic_section_template.Rmd" = get_default_ic_analysis_template(),
    "recovery_section_template.Rmd" = get_default_recovery_analysis_template(),
    "ppc_section_template.Rmd" = get_default_ppc_analysis_template(),
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
    "{{PPC_SECTION}}"
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
    "Parameter recovery analysis tests whether model parameters can be accurately recovered from simulated data.",
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
    "```"
  )
}

#' Default PPC analysis template
get_default_ppc_analysis_template <- function() {
  c(
    "# Posterior Predictive Checks Analysis",
    "",
    "Posterior predictive checks assess how well models reproduce observed behavioral patterns.",
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
