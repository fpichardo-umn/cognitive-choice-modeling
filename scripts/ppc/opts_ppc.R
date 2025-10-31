#!/usr/bin/env Rscript

#' Simple PPC Options - Just uncomment the opt list you want to test!
#' 
#' group parameter explanation:
#' - For pipeline scripts: "sing" or "hier" (fit type)
#' - For individual scripts: "batch_001" or "hier" (file naming)

# ========== CONFIGURATION ==========
opt_model = "rd_b1"
opt_task = "igt" 
opt_cohort = "adb"
opt_ses = "00"  # or "ses-01"
opt_group_name = "hier"  # for individual fits

# ========== PIPELINE SCRIPTS ==========

# Run full pipeline - Individual fits
# opt = list(
#   model = opt_model,
#   task = opt_task,
#   cohort = opt_cohort,
#   ses = opt_ses,
#   group = "sing",
#   group_name = opt_group_name,
#   n_sims = 25,
#   exclude_file = "Data/txt/subs/subject_ids_excludes.txt",
#   parallel = FALSE,
#   n_cores = 2,
#   sampling = "width",
#   width_control = 0.5,
#   rt_method = "remove",
#   RTbound_min_ms = 100,
#   RTbound_max_ms = 2500,
#   steps = "all",
#   force = FALSE,
#   render = TRUE,
#   ic_method = "loo"
# )

# Run full pipeline - Hierarchical fits
opt = list(
  model = opt_model,
  task = opt_task,
  cohort = opt_cohort,
  ses = opt_ses,
  group = "hier",
  group_name = opt_group_name,
  n_sims = 25,
  exclude_file = NULL,
  parallel = FALSE,
  n_cores = 2,
  sampling = "width",
  width_control = 0.5,
  rt_method = "mark",
  RTbound_min_ms = 50,
  RTbound_max_ms = 120000,
  steps = "all",
  force = FALSE,
  render = TRUE,
  ic_method = "loo"
)

# ========== INDIVIDUAL SCRIPTS ==========

# Run simulation - Individual fits
opt = list(
  model = opt_model,
  task = opt_task,
  cohort = opt_cohort,
  ses = opt_ses,
  group = "hier",
  group_name = opt_group_name,
  fit_file = NULL,
  n_sims = 25,
  exclude_file = NULL,#"Data/txt/subs/subject_ids_excludes.txt",
  output_dir = NULL,
  rt_method = "mark",
  RTbound_min_ms = 50,
  RTbound_max_ms = 120000,
  parallel = FALSE,
  n_cores = 2,
  sampling = "width",
  width_control = 0.5
)

# Run simulation - Hierarchical fits
opt = list(
  model = opt_model,
  task = opt_task,
  cohort = opt_cohort,
  ses = opt_ses,
  group = "hier",
  group_name = opt_group_name,
  fit_file = NULL,
  n_sims = 25,
  exclude_file = NULL, #"Data/txt/subs/subject_ids_excludes.txt",
  output_dir = NULL,
  rt_method = "remove",
  RTbound_min_ms = 100,
  RTbound_max_ms = 2500,
  parallel = FALSE,
  n_cores = 2,
  sampling = "width",
  width_control = 0.5
)

# Run stats - Individual fits
opt = list(
  model = opt_model,
  task = opt_task,
  cohort = opt_cohort,
  ses = opt_ses,
  group = opt_group_name,  # File naming: use batch name
  block_size = 20,
  output_dir = NULL,
  sim_file = NULL,
  exclude_file = "Data/txt/subs/subject_ids_excludes.txt"
)

# Run stats - Hierarchical fits
opt = list(
  model = opt_model,
  task = opt_task,
  cohort = opt_cohort,
  ses = opt_ses,
  group = "hier",  # File naming: use "hier"
  block_size = 20,
  output_dir = NULL,
  sim_file = NULL,
  exclude_file = "Data/txt/subs/subject_ids_excludes.txt"
)

# Run loglik - Individual fits
opt = list(
  model = opt_model,
  task = opt_task,
  cohort = opt_cohort,
  ses = opt_ses,
  group = opt_group_name,  # File naming: use batch name
  ic_method = "loo",
  output_dir = NULL,
  sim_file = NULL,
  rt_method = "mark",
  RTbound_min_ms = 50,
  RTbound_max_ms = 120000
)

# Run loglik - Hierarchical fits
opt = list(
  model = opt_model,
  task = opt_task,
  cohort = opt_cohort,
  ses = opt_ses,
  group = "hier",  # File naming: use "hier"
  ic_method = "loo",
  output_dir = NULL,
  sim_file = NULL,
  rt_method = "remove",
  RTbound_min_ms = 100,
  RTbound_max_ms = 2500
)

# Run report - Individual fits
opt = list(
  model = opt_model,
  task = opt_task,
  cohort = opt_cohort,
  ses = opt_ses,
  group = opt_group_name,  # File naming: use batch name
  output_dir = NULL,
  stats_file = NULL,
  force = FALSE,
  render = TRUE
)

# Run report - Hierarchical fits
opt = list(
  model = opt_model,
  task = opt_task,
  cohort = opt_cohort,
  ses = opt_ses,
  group = "hier",  # File naming: use "hier"
  output_dir = NULL,
  stats_file = NULL,
  force = FALSE,
  render = TRUE
)
