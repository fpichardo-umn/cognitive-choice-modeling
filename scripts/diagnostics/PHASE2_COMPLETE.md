# Phase 2: Analysis Modules - COMPLETE

## Status: ✓ COMPLETE

Date: 2025-01-09

## Components Implemented

### 1. Single Fit Diagnostics (`single_fit_diagnostics.R`)

**Main function:** `analyze_single_fit()`

**Sub-functions:**
- `extract_fit_metadata()` - Extract fit object metadata
- `analyze_convergence()` - Analyze R-hat, ESS, MCSE
- `analyze_sampling_diagnostics()` - Analyze divergences, treedepth, EBFMI
- `analyze_parameter_diagnostics()` - Per-parameter diagnostic table
- `prepare_single_fit_plots()` - Prepare plot data structures
- `get_problematic_parameters()` - Filter to problematic parameters only

**Returns:**
- metadata: Basic fit information
- convergence: R-hat, ESS, MCSE analysis
- sampling: Divergences, treedepth, EBFMI
- parameters: Per-parameter diagnostic table
- summary: Overall summary with status
- plot_data: Data for visualizations (optional)
- recommendations: Actionable recommendations

### 2. Batch Fit Diagnostics (`batch_fit_diagnostics.R`)

**Main function:** `analyze_batch_fits()`

**Sub-functions:**
- `aggregate_batch_diagnostics()` - Aggregate stats across subjects
- `aggregate_parameter_diagnostics()` - Parameter-level aggregation
- `identify_problematic_subjects_batch()` - Rank subjects by problems
- `determine_batch_status()` - Overall batch PASS/WARN/FAIL
- `create_batch_recommendations()` - Batch-specific recommendations
- `create_subject_table_for_export()` - Prepare CSV export

**Returns:**
- n_subjects: Number of subjects in batch
- subject_ids: List of subject IDs
- aggregate: Cross-subject statistics
- subjects: Individual subject summaries
- parameter_summary: Parameter-level aggregation
- problematic_subjects: Worst N subjects + all subjects table
- overall_status: Batch-level status
- recommendations: Actionable recommendations

### 3. Hierarchical Diagnostics (`hierarchical_diagnostics.R`)

**Main function:** `analyze_hierarchical_fit()`

**Sub-functions:**
- `identify_hierarchical_parameters()` - Parse parameter structure
- `analyze_subject_level_parameters()` - Subject parameter diagnostics
- `check_basic_hierarchical_health()` - Basic hyperparameter checks
- `check_shrinkage_basic()` - Simple shrinkage assessment (fit-based)
- `analyze_group_level_parameters()` - Full group parameter analysis (real hier)
- `analyze_shrinkage_full()` - Detailed shrinkage analysis (real hier)
- `determine_hierarchical_status()` - Overall status
- `create_hierarchical_recommendations()` - Hierarchical-specific recommendations

**Modes:**
- **Fit-based hierarchical** (default): Focus on subject estimates, basic checks
- **Real hierarchical** (--real_hier flag): Full group + subject analysis

**Returns:**
- fit_type: "fit_based_hierarchical" or "real_hierarchical"
- n_subjects: Number of subjects
- subject_analysis: Subject-level diagnostics
- basic_checks: Hyperparameter convergence, divergences
- shrinkage_check: Basic shrinkage assessment (fit-based)
- group_analysis: Group parameter diagnostics (real hier only)
- shrinkage_analysis: Full shrinkage analysis (real hier only)
- overall_status: Model-level status
- recommendations: Actionable recommendations

### 4. Problematic Subjects (`problematic_subjects.R`)

**Main function:** `identify_problematic_subjects()`

**Sub-functions:**
- `extract_subjects_from_batch()` - Extract from batch analysis
- `extract_subjects_from_hierarchical()` - Extract from hierarchical analysis
- `calculate_comprehensive_scores()` - Calculate/update problem scores
- `categorize_problems()` - Count by problem type
- `identify_specific_problem_types()` - Subjects by specific issues
- `create_cross_model_comparison()` - Compare across models
- `generate_problem_summary_text()` - Human-readable summary
- `export_problematic_subjects_csv()` - Export to CSV

**Returns:**
- summary: Overall statistics and problem breakdown
- all_subjects: Full subject table with scores
- worst_subjects: Top N worst subjects
- subjects_by_problem: Subjects grouped by problem type

## Key Features

### Comprehensive Diagnostics
- R-hat (convergence)
- ESS bulk and tail (effective sample size)
- MCSE (Monte Carlo standard error)
- Divergent transitions
- Max treedepth hits
- EBFMI (energy diagnostic)
- Energy distribution normality

### Flexible Analysis
- Works with single, batch, and hierarchical fits
- Automatic fit type detection
- Configurable thresholds
- Problem scoring with weights

### Hierarchical Support
- Two modes: fit-based vs real hierarchical
- Subject-level parameter focus
- Optional group-level analysis
- Shrinkage assessment

### Subject Ranking
- Weighted problem scores
- Multiple classification schemes
- Cross-model comparison support
- CSV export for large batches

## Integration

All analysis modules:
- Use `diagnostic_helpers.R` utilities
- Use `helper_diagnostics.R` existing functions
- Return structured lists (compatible with RDS)
- Accept threshold configuration
- Provide recommendations

## Testing

To test Phase 2 modules with actual fit data:

```r
library(here)
source(file.path(here::here(), "scripts", "diagnostics", "helpers", "diagnostic_helpers.R"))
source(file.path(here::here(), "scripts", "diagnostics", "analysis", "single_fit_diagnostics.R"))
source(file.path(here::here(), "scripts", "diagnostics", "analysis", "batch_fit_diagnostics.R"))

# Load a fit object
fit <- readRDS("path/to/fit/file.rds")

# Determine type
fit_type <- determine_fit_type(fit)
print(fit_type)

# Analyze based on type
if (fit_type == "single") {
  analysis <- analyze_single_fit(fit)
} else if (fit_type == "batch") {
  analysis <- analyze_batch_fits(fit)
}

print(analysis$overall_status)
print(analysis$recommendations)
```

## Next Steps

**Phase 3: Report Templates & Visualization**
- Create Rmd templates for each section
- Implement visualization functions
- Create educational content
- Build report generation system

Ready to proceed: **YES** ✓
