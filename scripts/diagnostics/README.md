# MCMC Diagnostics Pipeline

This directory contains scripts and configuration for automated MCMC diagnostic reporting.

## Structure

```
diagnostics/
├── run_diagnostics.R              # Main script (Phase 2+)
├── run_diagnostics_batch.R        # Batch processing (Phase 4+)
├── config/
│   ├── diagnostic_thresholds.yaml # Threshold definitions
│   └── report_config.yaml         # Report formatting options
├── helpers/
│   ├── diagnostic_dirs.R          # Directory management [COMPLETE]
│   └── diagnostic_helpers.R       # Utility functions [COMPLETE]
├── analysis/                      # Phase 2
│   ├── single_fit_diagnostics.R
│   ├── batch_fit_diagnostics.R
│   ├── hierarchical_diagnostics.R
│   └── problematic_subjects.R
├── reports/                       # Phase 3
│   ├── generate_diagnostic_report.R
│   └── templates/
│       ├── main_diagnostic_template.Rmd
│       ├── single_fit_section.Rmd
│       ├── batch_fit_section.Rmd
│       ├── hierarchical_section.Rmd
│       ├── problematic_subjects_section.Rmd
│       └── diagnostic_explanations.Rmd
└── visualization/                 # Phase 3
    ├── convergence_plots.R
    ├── diagnostic_summary_plots.R
    └── subject_comparison_plots.R
```

## Phase 1: Core Infrastructure [COMPLETE]

### Completed Components

1. **Helper Functions** (`helpers/`)
   - `diagnostic_dirs.R`: Directory management and path resolution
   - `diagnostic_helpers.R`: Diagnostic processing utilities

2. **Configuration Files** (`config/`)
   - `diagnostic_thresholds.yaml`: Thresholds for PASS/WARN/FAIL classification
   - `report_config.yaml`: Report formatting and display options

### Key Functions

#### diagnostic_dirs.R

```r
# Get output directories
get_diagnostics_base_dir(task, cohort, session)
get_diagnostics_reports_dir(task, cohort, session)
get_diagnostics_summaries_dir(task, cohort, session)
get_diagnostics_tables_dir(task, cohort, session)
get_diagnostics_logs_dir(task, cohort, session)

# Ensure directories exist
ensure_diagnostics_dirs(task, cohort, session)

# Get file paths with BIDS naming
get_diagnostics_report_path(task, cohort, session, group, model)
get_diagnostics_summary_path(task, cohort, session, group, model)
get_diagnostics_table_path(task, cohort, session, group, model)

# Find fit files automatically
find_fit_file(task, cohort, session, group, model, subid)
```

#### diagnostic_helpers.R

```r
# Load configurations
load_diagnostic_thresholds()
load_report_config()

# Classify diagnostics
classify_diagnostic(value, thresholds, metric_type)
aggregate_status(classifications)

# Problem scoring
calculate_problem_score(subject_diagnostics, thresholds)

# Fit type detection
determine_fit_type(fit)
extract_subject_ids(fit)

# Summary creation
create_single_fit_summary(fit, subject_id, thresholds)

# Formatting
format_status_badge(status)
create_recommendations(diagnostic_summary, fit_type)
```

### Configuration

#### Diagnostic Thresholds

- **R-hat**: excellent ≤ 1.01, acceptable ≤ 1.05, problematic > 1.1
- **ESS ratio**: good ≥ 0.7, acceptable ≥ 0.5, problematic < 0.3
- **Divergence rate**: acceptable ≤ 0.1%, borderline ≤ 1%, problematic > 5%
- **MCSE ratio**: acceptable ≤ 5%, borderline ≤ 10%, problematic > 15%
- **EBFMI**: acceptable ≥ 0.3, problematic < 0.2

#### Problem Weights

- High R-hat: 5
- Divergences: 3
- Low ESS: 2
- High MCSE: 1
- Low EBFMI: 2

### Output Location

```
Outputs/{TASK}/fits/diagnostics/{COHORT}/[ses-{SESSION}/]
├── reports/      # HTML reports
├── summaries/    # RDS diagnostic summaries
├── tables/       # CSV tables (for large batches)
└── logs/         # Processing logs
```

### File Naming (BIDS-style)

**HTML Report:**
```
task-{task}_cohort-{cohort}_[ses-{ses}_]group-{group}_model-{model}_desc-diagnostics.html
```

**RDS Summary:**
```
task-{task}_cohort-{cohort}_[ses-{ses}_]group-{group}_model-{model}_desc-diagsummary.rds
```

**CSV Table:**
```
task-{task}_cohort-{cohort}_[ses-{ses}_]group-{group}_model-{model}_desc-allsubjects.csv
```

## Testing Phase 1

To test the completed Phase 1 components:

```r
# Test directory creation
library(here)
source(file.path(here::here(), "scripts", "diagnostics", "helpers", "diagnostic_dirs.R"))

dirs <- ensure_diagnostics_dirs("igt_mod", "adb", "01")
print(dirs)

# Test configuration loading
source(file.path(here::here(), "scripts", "diagnostics", "helpers", "diagnostic_helpers.R"))

thresholds <- load_diagnostic_thresholds()
config <- load_report_config()

print(thresholds$thresholds$rhat)
print(config$report$n_worst_subjects)

# Test classification
rhat_value <- 1.03
status <- classify_diagnostic(rhat_value, thresholds, "rhat")
print(paste("R-hat", rhat_value, "is:", status))

# Test status badge formatting
print(format_status_badge("PASS"))
print(format_status_badge("WARN"))
print(format_status_badge("FAIL"))
```

## Next Steps

### Phase 2: Analysis Modules
- Implement single fit diagnostics extraction
- Implement batch fit aggregation
- Implement hierarchical fit analysis
- Implement problematic subject identification

### Phase 3: Report Templates
- Create Rmd templates for each section
- Create educational content
- Create visualization functions

### Phase 4: Orchestration
- Implement main run_diagnostics.R script
- Implement batch processing script
- Add comprehensive error handling
- Add logging

### Phase 5: Testing & Documentation
- Test with actual fit objects
- Generate example reports
- Create user documentation
- Add error handling edge cases

## Dependencies

- R packages: `here`, `yaml`, `dplyr`, `posterior`, `ggplot2`, `bayesplot`, `knitr`, `rmarkdown`
- Existing scripts: `helper_diagnostics.R`, `helper_dirs.R`, `helper_common.R`, `helper_functions_cmdSR.R`

## Notes

- All paths use BIDS-style naming conventions
- Configurations are YAML for easy editing
- Helper functions are modular and reusable
- Design follows model_comparison pipeline pattern
