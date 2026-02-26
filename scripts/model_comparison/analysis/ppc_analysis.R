#!/usr/bin/env Rscript

#' PPC Analysis by Behavioral Domains
#' @description Analyze posterior predictive check performance organized by behavioral domains

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(here)
})

# Load helper functions
source(file.path(here::here(), "scripts", "model_comparison", "helpers", "model_comparison_helpers.R"))
source(file.path(here::here(), "scripts", "ppc", "helpers", "task_config.R"))

#' Analyze PPC performance for all models organized by behavioral domains
#' @param comparison_data List of model data from load_comparison_data
#' @param models_by_type List organizing models by type  
#' @param task Task name for task-specific analysis
#' @return List with PPC analysis results
analyze_ppc_by_groups <- function(comparison_data, models_by_type, task) {
  message("Analyzing PPC performance by behavioral domains...")
  
  # Initialize results
  results <- list(
    by_domain_and_model = data.frame(),
    domain_summary = data.frame(),
    model_summary = data.frame(),
    extreme_failures = data.frame(),
    behavioral_patterns = list(),
    task_specific_analysis = list()
  )
  
  # Get task configuration for behavioral domain definitions
  task_config <- get_task_config(task)
  behavioral_domains <- define_behavioral_domains(task)
  
  # Process each model
  all_domain_model_data <- list()
  
  for (model_name in names(comparison_data)) {
    model_data <- comparison_data[[model_name]]
    
    if (is.null(model_data$ppc)) {
      warning("No PPC data for model: ", model_name)
      next
    }
    
    # Analyze this model's PPC by domains
    model_analysis <- analyze_single_model_ppc(model_data, model_name, behavioral_domains, task)
    
    if (!is.null(model_analysis)) {
      all_domain_model_data[[model_name]] <- model_analysis
    }
  }
  
  if (length(all_domain_model_data) == 0) {
    warning("No PPC data available for any models")
    return(results)
  }
  
  # Combine results across models
  results$by_domain_and_model <- do.call(rbind, lapply(names(all_domain_model_data), function(model) {
    data <- all_domain_model_data[[model]]$by_domain
    data$model <- model
    data$model_type <- classify_model_type(model, task)
    return(data)
  }))
  
  # Create domain-level summary (across models)
  if (nrow(results$by_domain_and_model) > 0) {
    results$domain_summary <- results$by_domain_and_model %>%
      group_by(domain) %>%
      summarise(
        n_models = n_distinct(model),
        mean_ppp = mean(mean_ppp, na.rm = TRUE),
        median_ppp = median(mean_ppp, na.rm = TRUE),
        proportion_extreme = mean(proportion_extreme, na.rm = TRUE),
        n_statistics = sum(n_statistics, na.rm = TRUE),
        mean_deviation = mean(mean_deviation, na.rm = TRUE),
        domain_quality = classify_ppc_quality(mean(mean_ppp, na.rm = TRUE), mean(proportion_extreme, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      arrange(proportion_extreme, abs(mean_ppp - 0.5))
    
    # ──────────────────────────────────────────────────────────────────────
    # FIX: Create model-level summary using PROPERLY WEIGHTED metrics.
    #
    # Previously this used mean(proportion_extreme) across domains, which
    # gave equal weight to small and large domains — a domain with 50
    # statistics and 30% extreme rate would count the same as one with
    # 1000 statistics and 1% extreme rate, massively inflating the result.
    #
    # Now we use sum(n_extreme) / sum(n_statistics) for the true overall
    # rate, and weighted.mean() for ppp and deviation, matching what the
    # individual model PPC reports compute (total extreme / total stats).
    # ──────────────────────────────────────────────────────────────────────
    results$model_summary <- results$by_domain_and_model %>%
      group_by(model, model_type) %>%
      summarise(
        n_domains = n_distinct(domain),
        overall_mean_ppp = weighted.mean(mean_ppp, n_statistics, na.rm = TRUE),
        overall_proportion_extreme = sum(n_extreme, na.rm = TRUE) / sum(n_statistics, na.rm = TRUE),
        total_statistics = sum(n_statistics, na.rm = TRUE),
        overall_deviation = weighted.mean(mean_deviation, n_statistics, na.rm = TRUE),
        model_quality = classify_ppc_quality(
          weighted.mean(mean_ppp, n_statistics, na.rm = TRUE),
          sum(n_extreme, na.rm = TRUE) / sum(n_statistics, na.rm = TRUE)
        ),
        worst_domain = domain[which.max(proportion_extreme)][1],
        best_domain = domain[which.min(proportion_extreme)][1],
        .groups = "drop"
      ) %>%
      arrange(overall_proportion_extreme, abs(overall_mean_ppp - 0.5))
  }
  
  # Identify extreme failures across models
  results$extreme_failures <- identify_extreme_ppc_failures(all_domain_model_data, behavioral_domains, task)
  
  # Analyze behavioral patterns
  results$behavioral_patterns <- analyze_behavioral_patterns(all_domain_model_data, models_by_type, task)
  
  # Task-specific analysis
  results$task_specific_analysis <- conduct_task_specific_ppc_analysis(all_domain_model_data, task, task_config)
  
  # Add metadata
  results$metadata <- list(
    n_models_analyzed = length(all_domain_model_data),
    models_analyzed = names(all_domain_model_data),
    domains_analyzed = names(behavioral_domains),
    task = task,
    analysis_timestamp = Sys.time()
  )
  
  message("PPC analysis complete for ", length(all_domain_model_data), " models")
  return(results)
}

#' Define behavioral domains for a task
#' @param task Task name
#' @return List of behavioral domain definitions
define_behavioral_domains <- function(task) {
  
  if (task == "igt_mod") {
    # Modified IGT behavioral domains
    return(list(
      choice_patterns = list(
        description = "Overall choice behavior and deck preferences",
        statistics = c("play_ratio", "pass_ratio", "play_ratio_deck1", "play_ratio_deck2", 
                      "play_ratio_deck3", "play_ratio_deck4"),
        session_types = "session"
      ),
      performance = list(
        description = "Learning and performance metrics",
        statistics = c("good_play_ratio", "bad_play_ratio", "net_score", "mean_earnings", "total_earnings"),
        session_types = c("session", "block")
      ),
      learning_dynamics = list(
        description = "Learning curves and temporal patterns",
        statistics = c("good_play_ratio", "bad_play_ratio", "net_score"),
        session_types = "block"
      ),
      rt_patterns = list(
        description = "Response time distributions and patterns",
        statistics = c("rt_mean", "rt_sd", "rt_q10", "rt_q50", "rt_q90", 
                      "rt_mean_play", "rt_mean_pass", "rt_mean_good", "rt_mean_bad"),
        session_types = c("session", "block")
      )
    )
  )
  } else if (task == "igt") {
    # Original IGT behavioral domains
    return(list(
      choice_patterns = list(
        description = "Deck selection patterns and preferences",
        statistics = c("deck1_freq", "deck2_freq", "deck3_freq", "deck4_freq"),
        session_types = "session"
      ),
      performance = list(
        description = "Learning performance and strategy",
        statistics = c("good_deck_freq", "bad_deck_freq", "net_score", "mean_earnings", "total_earnings"),
        session_types = c("session", "block")
      ),
      strategies = list(
        description = "Decision strategies and perseveration",
        statistics = c("win_stay", "lose_shift", "perseveration"),
        session_types = c("session", "block")
      ),
      learning_dynamics = list(
        description = "Learning curves and temporal patterns", 
        statistics = c("good_deck_freq", "bad_deck_freq", "net_score"),
        session_types = "block"
      ),
      rt_patterns = list(
        description = "Response time distributions and patterns",
        statistics = c("rt_mean", "rt_sd", "rt_q10", "rt_q50", "rt_q90",
                      "rt_deck1", "rt_deck2", "rt_deck3", "rt_deck4"),
        session_types = c("session", "block")
      )
    )
  )
  } else {
    stop("Unknown task for behavioral domain definition: ", task)
  }
}

#' Analyze PPC performance for a single model
#' @param model_data Model data including PPC results
#' @param model_name Name of the model
#' @param behavioral_domains Behavioral domain definitions
#' @param task Task name
#' @return List with PPC analysis for this model
analyze_single_model_ppc <- function(model_data, model_name, behavioral_domains, task) {
  ppc_data <- model_data$ppc
  
  if (is.null(ppc_data) || nrow(ppc_data) == 0) {
    return(NULL)
  }
  
  # Calculate domain-level metrics
  domain_metrics <- list()
  
  for (domain_name in names(behavioral_domains)) {
    domain_def <- behavioral_domains[[domain_name]]
    
    # Filter PPC data for this domain
    domain_data <- ppc_data %>%
      filter(
        statistic %in% domain_def$statistics,
        (session %in% domain_def$session_types | 
         grepl(paste(domain_def$session_types, collapse = "|"), session))
      )
    
    if (nrow(domain_data) == 0) next
    
    # Calculate domain-level metrics
    domain_metrics[[domain_name]] <- calculate_domain_ppc_metrics(domain_data, domain_name, model_name)
  }
  
  # Combine domain metrics
  by_domain <- do.call(rbind, domain_metrics)
  
  # Calculate overall model PPC metrics
  overall_metrics <- calculate_overall_ppc_metrics(ppc_data, model_name)
  
  return(list(
    by_domain = by_domain,
    overall = overall_metrics,
    raw_data = ppc_data
  ))
}

#' Calculate PPC metrics for a behavioral domain
#' @param domain_data PPC data filtered for this domain
#' @param domain_name Name of the behavioral domain
#' @param model_name Name of the model
#' @return Data frame with domain-level PPC metrics
calculate_domain_ppc_metrics <- function(domain_data, domain_name, model_name) {
  if (nrow(domain_data) == 0) {
    return(data.frame())
  }
  
  # Calculate summary metrics
  mean_ppp <- mean(domain_data$ppp, na.rm = TRUE)
  median_ppp <- median(domain_data$ppp, na.rm = TRUE)
  sd_ppp <- sd(domain_data$ppp, na.rm = TRUE)
  
  # Count extreme values
  n_extreme <- sum(domain_data$extreme_ppp %in% TRUE, na.rm = TRUE)
  proportion_extreme <- n_extreme / nrow(domain_data)
  
  # Calculate deviation from ideal (PPP = 0.5)
  mean_deviation <- mean(abs(domain_data$ppp - 0.5), na.rm = TRUE)
  
  # Count statistics
  n_statistics <- nrow(domain_data)
  n_subjects <- n_distinct(domain_data$subject_id)
  
  # Identify worst performing statistics
  worst_stats <- domain_data %>%
    arrange(desc(extreme_ppp), desc(abs(ppp - 0.5))) %>%
    head(3) %>%
    pull(statistic)
  
  return(data.frame(
    domain = domain_name,
    mean_ppp = mean_ppp,
    median_ppp = median_ppp,
    sd_ppp = sd_ppp,
    n_extreme = n_extreme,
    proportion_extreme = proportion_extreme,
    mean_deviation = mean_deviation,
    n_statistics = n_statistics,
    n_subjects = n_subjects,
    worst_statistics = paste(worst_stats, collapse = ", "),
    domain_quality = classify_ppc_quality(mean_ppp, proportion_extreme),
    stringsAsFactors = FALSE
  ))
}

#' Calculate overall PPC metrics for a model
#' @param ppc_data Full PPC data for the model
#' @param model_name Name of the model
#' @return Data frame with overall PPC metrics
calculate_overall_ppc_metrics <- function(ppc_data, model_name) {
  if (nrow(ppc_data) == 0) {
    return(data.frame())
  }
  
  # Overall metrics
  overall_mean_ppp <- mean(ppc_data$ppp, na.rm = TRUE)
  overall_proportion_extreme <- mean(ppc_data$extreme_ppp %in% TRUE, na.rm = TRUE)
  overall_deviation <- mean(abs(ppc_data$ppp - 0.5), na.rm = TRUE)
  
  # By category
  by_category <- ppc_data %>%
    group_by(category) %>%
    summarise(
      mean_ppp = mean(ppp, na.rm = TRUE),
      proportion_extreme = mean(extreme_ppp %in% TRUE, na.rm = TRUE),
      n_statistics = n(),
      .groups = "drop"
    )
  
  return(list(
    overall = data.frame(
      overall_mean_ppp = overall_mean_ppp,
      overall_proportion_extreme = overall_proportion_extreme,
      overall_deviation = overall_deviation,
      total_statistics = nrow(ppc_data),
      model_quality = classify_ppc_quality(overall_mean_ppp, overall_proportion_extreme)
    ),
    by_category = by_category
  ))
}

#' Classify PPC quality based on PPP values and extreme proportion
#' @param mean_ppp Mean PPP value
#' @param proportion_extreme Proportion of extreme PPP values
#' @return Quality label
classify_ppc_quality <- function(mean_ppp, proportion_extreme) {
  if (is.na(mean_ppp) || is.na(proportion_extreme)) return("unknown")
  
  deviation_from_ideal <- abs(mean_ppp - 0.5)
  
  if (proportion_extreme > 0.2) return("poor")
  if (proportion_extreme > 0.1 && deviation_from_ideal > 0.2) return("concerning")
  if (proportion_extreme > 0.05 || deviation_from_ideal > 0.15) return("acceptable")
  if (deviation_from_ideal < 0.1 && proportion_extreme < 0.05) return("excellent")
  return("good")
}

#' Identify extreme PPC failures across models
#' @param all_domain_model_data All model PPC analyses
#' @param behavioral_domains Domain definitions
#' @param task Task name (required for classify_model_type)
#' @return Data frame with extreme failures
identify_extreme_ppc_failures <- function(all_domain_model_data, behavioral_domains, task) {
  extreme_failures <- list()
  
  for (model_name in names(all_domain_model_data)) {
    model_data <- all_domain_model_data[[model_name]]
    
    if (is.null(model_data$raw_data)) next
    
    # Find extreme failures in raw data
    extreme_data <- model_data$raw_data %>%
      filter(extreme_ppp %in% TRUE) %>%
      mutate(
        model = model_name,
        model_type = classify_model_type(model_name, task),
        failure_severity = case_when(
          ppp < 0.01 | ppp > 0.99 ~ "severe",
          ppp < 0.025 | ppp > 0.975 ~ "moderate",
          TRUE ~ "mild"
        )
      )
    
    if (nrow(extreme_data) > 0) {
      extreme_failures[[model_name]] <- extreme_data
    }
  }
  
  if (length(extreme_failures) == 0) {
    return(data.frame())
  }
  
  # Combine and summarize
  combined_failures <- do.call(rbind, extreme_failures)
  
  return(combined_failures)
}

#' Analyze behavioral patterns across model types
#' @param all_domain_model_data All model PPC analyses
#' @param models_by_type Models organized by type
#' @param task Task name
#' @return List with behavioral pattern analysis
analyze_behavioral_patterns <- function(all_domain_model_data, models_by_type, task) {
  patterns <- list()
  
  # Combine all domain data
  all_domain_data <- do.call(rbind, lapply(names(all_domain_model_data), function(model) {
    data <- all_domain_model_data[[model]]$by_domain
    if (!is.null(data) && nrow(data) > 0) {
      data$model <- model
      data$model_type <- classify_model_type(model, task)
      return(data)
    }
    return(data.frame())
  }))
  
  if (nrow(all_domain_data) == 0) {
    return(patterns)
  }
  
  # Pattern 1: Which domains are hardest to capture across model types?
  patterns$difficult_domains <- all_domain_data %>%
    group_by(domain, model_type) %>%
    summarise(
      n_models = n_distinct(model),
      mean_proportion_extreme = mean(proportion_extreme, na.rm = TRUE),
      mean_deviation = mean(mean_deviation, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_proportion_extreme))
  
  # Pattern 2: Do hybrid models perform better on RT patterns?
  if ("rt_patterns" %in% all_domain_data$domain) {
    patterns$rt_performance_by_type <- all_domain_data %>%
      filter(domain == "rt_patterns") %>%
      group_by(model_type) %>%
      summarise(
        n_models = n_distinct(model),
        mean_proportion_extreme = mean(proportion_extreme, na.rm = TRUE),
        mean_deviation = mean(mean_deviation, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # Pattern 3: Learning dynamics performance
  if ("learning_dynamics" %in% all_domain_data$domain) {
    patterns$learning_performance_by_type <- all_domain_data %>%
      filter(domain == "learning_dynamics") %>%
      group_by(model_type) %>%
      summarise(
        n_models = n_distinct(model),
        mean_proportion_extreme = mean(proportion_extreme, na.rm = TRUE),
        mean_deviation = mean(mean_deviation, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # Pattern 4: Best and worst performing models by domain
  patterns$model_performance_by_domain <- all_domain_data %>%
    group_by(domain) %>%
    arrange(proportion_extreme, mean_deviation) %>%
    summarise(
      best_model = first(model),
      best_model_extreme_prop = first(proportion_extreme),
      worst_model = last(model), 
      worst_model_extreme_prop = last(proportion_extreme),
      performance_range = last(proportion_extreme) - first(proportion_extreme),
      .groups = "drop"
    )
  
  return(patterns)
}

#' Conduct task-specific PPC analysis
#' @param all_domain_model_data All model PPC analyses
#' @param task Task name
#' @param task_config Task configuration
#' @return List with task-specific analysis
conduct_task_specific_ppc_analysis <- function(all_domain_model_data, task, task_config) {
  task_analysis <- list()
  
  if (task == "igt_mod") {
    # mIGT-specific analysis
    task_analysis$play_pass_analysis <- analyze_play_pass_patterns(all_domain_model_data)
    task_analysis$deck_preference_analysis <- analyze_deck_preference_patterns(all_domain_model_data)
  } else if (task == "igt") {
    # IGT-specific analysis
    task_analysis$deck_selection_analysis <- analyze_deck_selection_patterns(all_domain_model_data)
    task_analysis$strategy_analysis <- analyze_strategy_patterns(all_domain_model_data)
  }
  
  return(task_analysis)
}

#' Analyze play/pass patterns for mIGT
#' @param all_domain_model_data All model PPC analyses  
#' @return Data frame with play/pass analysis
analyze_play_pass_patterns <- function(all_domain_model_data) {
  # Extract play/pass related statistics
  play_pass_stats <- c("play_ratio", "pass_ratio", "good_play_ratio", "bad_play_ratio")
  
  results <- list()
  
  for (model_name in names(all_domain_model_data)) {
    model_data <- all_domain_model_data[[model_name]]$raw_data
    
    if (is.null(model_data)) next
    
    play_pass_data <- model_data %>%
      filter(statistic %in% play_pass_stats) %>%
      group_by(statistic) %>%
      summarise(
        mean_ppp = mean(ppp, na.rm = TRUE),
        proportion_extreme = mean(extreme_ppp %in% TRUE, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(model = model_name)
    
    if (nrow(play_pass_data) > 0) {
      results[[model_name]] <- play_pass_data
    }
  }
  
  if (length(results) == 0) return(data.frame())
  
  return(do.call(rbind, results))
}

#' Analyze deck preference patterns
#' @param all_domain_model_data All model PPC analyses
#' @return Data frame with deck preference analysis
analyze_deck_preference_patterns <- function(all_domain_model_data) {
  deck_stats <- c("play_ratio_deck1", "play_ratio_deck2", "play_ratio_deck3", "play_ratio_deck4")
  
  results <- list()
  
  for (model_name in names(all_domain_model_data)) {
    model_data <- all_domain_model_data[[model_name]]$raw_data
    
    if (is.null(model_data)) next
    
    deck_data <- model_data %>%
      filter(statistic %in% deck_stats) %>%
      group_by(statistic) %>%
      summarise(
        mean_ppp = mean(ppp, na.rm = TRUE),
        proportion_extreme = mean(extreme_ppp %in% TRUE, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        model = model_name,
        deck = as.numeric(gsub("play_ratio_deck", "", statistic)),
        deck_type = ifelse(deck %in% c(3, 4), "good", "bad")
      )
    
    if (nrow(deck_data) > 0) {
      results[[model_name]] <- deck_data
    }
  }
  
  if (length(results) == 0) return(data.frame())
  
  return(do.call(rbind, results))
}

#' Analyze deck selection patterns for IGT
#' @param all_domain_model_data All model PPC analyses
#' @return Data frame with deck selection analysis  
analyze_deck_selection_patterns <- function(all_domain_model_data) {
  deck_stats <- c("deck1_freq", "deck2_freq", "deck3_freq", "deck4_freq", "good_deck_freq", "bad_deck_freq")
  
  results <- list()
  
  for (model_name in names(all_domain_model_data)) {
    model_data <- all_domain_model_data[[model_name]]$raw_data
    
    if (is.null(model_data)) next
    
    deck_data <- model_data %>%
      filter(statistic %in% deck_stats) %>%
      group_by(statistic) %>%
      summarise(
        mean_ppp = mean(ppp, na.rm = TRUE),
        proportion_extreme = mean(extreme_ppp %in% TRUE, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(model = model_name)
    
    if (nrow(deck_data) > 0) {
      results[[model_name]] <- deck_data
    }
  }
  
  if (length(results) == 0) return(data.frame())
  
  return(do.call(rbind, results))
}

#' Analyze strategy patterns for IGT
#' @param all_domain_model_data All model PPC analyses
#' @return Data frame with strategy analysis
analyze_strategy_patterns <- function(all_domain_model_data) {
  strategy_stats <- c("win_stay", "lose_shift", "perseveration")
  
  results <- list()
  
  for (model_name in names(all_domain_model_data)) {
    model_data <- all_domain_model_data[[model_name]]$raw_data
    
    if (is.null(model_data)) next
    
    strategy_data <- model_data %>%
      filter(statistic %in% strategy_stats) %>%
      group_by(statistic) %>%
      summarise(
        mean_ppp = mean(ppp, na.rm = TRUE),
        proportion_extreme = mean(extreme_ppp %in% TRUE, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(model = model_name)
    
    if (nrow(strategy_data) > 0) {
      results[[model_name]] <- strategy_data
    }
  }
  
  if (length(results) == 0) return(data.frame())
  
  return(do.call(rbind, results))
}
