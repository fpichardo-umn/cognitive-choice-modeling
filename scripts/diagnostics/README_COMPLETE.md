# MCMC Diagnostics Pipeline - Complete Documentation

## Overview

A comprehensive pipeline for automated MCMC diagnostic analysis and reporting. Analyzes single subjects, batches, and hierarchical models with publication-ready HTML reports.

**Status:** ✓ Production Ready

---

## Quick Start

### Basic Usage

```bash
# Single subject
Rscript scripts/diagnostics/run_diagnostics.R \
  --task igt_mod --cohort adb --model ev_ddm --group sing --subid 0001

# Batch of subjects
Rscript scripts/diagnostics/run_diagnostics.R \
  --task igt_mod --cohort adb --model ev_ddm --group batch_001

# Hierarchical model
Rscript scripts/diagnostics/run_diagnostics.R \
  --task igt_mod --cohort adb --model ev_ddm --group hier

# Multiple models at once
Rscript scripts/diagnostics/run_diagnostics_batch.R \
  --task igt_mod --cohort adb --models all --parallel
```

---

## Features

### Comprehensive Diagnostics
- **R-hat** - Convergence assessment (Gelman-Rubin statistic)
- **ESS** - Effective sample size (bulk and tail)
- **MCSE** - Monte Carlo standard error
- **Divergences** - Problematic transitions detection
- **Tree Depth** - NUTS-specific diagnostics
- **EBFMI** - Energy distribution assessment

### Three Fit Types Supported
1. **Single Subject** - Detailed individual analysis
2. **Batch** - Aggregate analysis across multiple subjects
3. **Hierarchical** - Subject + group level analysis (two modes)

### Intelligent Reporting
- **Status Classification** - PASS/WARN/FAIL badges
- **Problem Scoring** - Weighted severity assessment
- **Actionable Recommendations** - What to do about issues
- **Educational Content** - Explains each diagnostic metric
- **Publication-Ready Plots** - Professional visualizations

### Flexible Configuration
- **Configurable Thresholds** - YAML-based customization
- **Multiple Output Formats** - HTML, RDS, CSV
- **Batch Processing** - Sequential or parallel
- **Auto-Detection** - Fit files and types

---

## Installation

### Dependencies

R packages (automatically loaded):
```r
# Core
optparse, here

# Analysis
cmdstanr, posterior, bayesplot

# Visualization
ggplot2, dplyr, tidyr, gridExtra

# Reporting
rmarkdown, knitr
```

### Setup

1. Ensure scripts/diagnostics/ exists in your project
2. Verify helper scripts are available:
   - `scripts/helpers/helper_functions_cmdSR.R`
   - `scripts/helpers/helper_common.R`
   - `scripts/helpers/helper_dirs.R`
   - `scripts/helpers/helper_diagnostics.R`

3. Create output directories (automatic on first run)

---

## Usage

### Command Line Arguments

#### `run_diagnostics.R`

**Required:**
- `--task/-k` - Task name (e.g., igt_mod)
- `--cohort/-c` - Cohort identifier (e.g., adb)
- `--model/-m` - Model name (e.g., ev_ddm)
- `--group/-g` - Group type (sing, batch_XXX, hier)

**For Single Fits:**
- `--subid` - Subject ID (required if group=sing)

**Optional:**
- `--session/-s` - Session identifier
- `--fit_file/-f` - Explicit fit file path (overrides auto-detection)
- `--output_dir/-o` - Custom output directory
- `--render_html` - Generate HTML report (default: TRUE)
- `--save_rds` - Save RDS summary (default: TRUE)
- `--table_threshold` - Max subjects in HTML table (default: 50)
- `--n_worst_subs` - Number of worst subjects to highlight (default: 10)
- `--real_hier` - Full hierarchical analysis (default: FALSE)
- `--force` - Overwrite existing outputs
- `--verbose/-v` - Detailed progress messages
- `--dry_run` - Preview without executing

#### `run_diagnostics_batch.R`

**Required:**
- `--task/-k` - Task name
- `--cohort/-c` - Cohort identifier

**Optional:**
- `--session/-s` - Session identifier
- `--models/-m` - 'all' or comma-separated list (default: all)
- `--groups/-g` - 'all' or comma-separated list (default: all)
- `--parallel/-p` - Enable parallel processing
- `--n_cores` - Number of cores (default: 4)
- Plus all options from run_diagnostics.R

---

## Examples

### 1. Single Subject Analysis

```bash
Rscript scripts/diagnostics/run_diagnostics.R \
  --task igt_mod \
  --cohort adb \
  --session 01 \
  --model ev_ddm \
  --group sing \
  --subid 0042 \
  --verbose
```

**Output:**
- HTML report with detailed diagnostics
- Per-parameter convergence metrics
- Visual diagnostics (traceplots, densities)
- Recommendations

### 2. Batch Analysis

```bash
Rscript scripts/diagnostics/run_diagnostics.R \
  --task igt_mod \
  --cohort adb \
  --model pvl \
  --group batch_001
```

**Output:**
- HTML report with aggregate statistics
- Distribution of diagnostics across subjects
- Top 10 worst subjects table
- CSV file with all subjects (if N > 50)
- Parameter-level problem analysis

### 3. Hierarchical Model (Fit-Based)

```bash
Rscript scripts/diagnostics/run_diagnostics.R \
  --task igt_mod \
  --cohort adb \
  --model ev_ddm \
  --group hier
```

**Output:**
- Subject-level diagnostics (focus)
- Basic hyperparameter checks
- Shrinkage assessment
- Worst subjects highlighted

### 4. Hierarchical Model (Real)

```bash
Rscript scripts/diagnostics/run_diagnostics.R \
  --task igt_mod \
  --cohort adb \
  --model ev_ddm \
  --group hier \
  --real_hier
```

**Output:**
- Full subject-level diagnostics
- Complete group-level analysis
- Detailed shrinkage analysis
- Hyperparameter convergence tracking

### 5. Batch Processing

```bash
# All models
Rscript scripts/diagnostics/run_diagnostics_batch.R \
  --task igt_mod \
  --cohort adb \
  --models all \
  --groups all \
  --parallel \
  --n_cores 8

# Specific models
Rscript scripts/diagnostics/run_diagnostics_batch.R \
  --task igt_mod \
  --cohort adb \
  --models ev,pvl,ddm \
  --groups batch_001,batch_002
```

**Output:**
- Individual reports for each model/group
- Index page linking all reports
- Batch processing log
- Success/failure summary

### 6. Custom Configuration

```bash
Rscript scripts/diagnostics/run_diagnostics.R \
  --task igt_mod \
  --cohort adb \
  --model ev_ddm \
  --group batch_001 \
  --table_threshold 100 \
  --n_worst_subs 20
```

### 7. Dry Run

```bash
Rscript scripts/diagnostics/run_diagnostics.R \
  --task igt_mod \
  --cohort adb \
  --model ev_ddm \
  --group batch_001 \
  --dry_run
```

---

## Output Structure

```
Outputs/{TASK}/fits/diagnostics/{COHORT}/[ses-{SESSION}/]
├── reports/
│   ├── task-igt_mod_cohort-adb_group-batch_001_model-ev_ddm_desc-diagnostics.html
│   ├── task-igt_mod_cohort-adb_group-hier_model-pvl_desc-diagnostics.html
│   └── index.html  # Generated by batch processing
│
├── summaries/
│   ├── task-igt_mod_cohort-adb_group-batch_001_model-ev_ddm_desc-diagsummary.rds
│   └── task-igt_mod_cohort-adb_group-hier_model-pvl_desc-diagsummary.rds
│
├── tables/
│   └── task-igt_mod_cohort-adb_group-batch_001_model-ev_ddm_desc-allsubjects.csv
│
└── logs/
    ├── task-igt_mod_cohort-adb_group-batch_001_model-ev_ddm_desc-diagnostics_*.log
    └── batch_diagnostics_*.log
```

---

## Configuration

### Threshold Configuration (`config/diagnostic_thresholds.yaml`)

```yaml
thresholds:
  rhat:
    excellent: 1.01
    acceptable: 1.05
    problematic: 1.1
  ess_ratio:
    good: 0.7
    acceptable: 0.5
    problematic: 0.3
  # ... etc
```

### Report Configuration (`config/report_config.yaml`)

```yaml
report:
  n_worst_subjects: 10
  table_threshold: 50
  plots:
    traceplots:
      n_params_sample: 10
  # ... etc
```

Edit these files to customize behavior.

---

## Report Contents

### HTML Report Sections

1. **Executive Summary**
   - Overall status badge (PASS/WARN/FAIL)
   - Key findings
   - Quick assessment
   - Recommendations

2. **Understanding MCMC Diagnostics**
   - Comprehensive explanations
   - What each metric means
   - What to look for
   - What to do about problems
   - Links to Stan documentation

3. **Diagnostic Results**
   - Fit-type-specific analysis
   - Convergence metrics (R-hat, ESS, MCSE)
   - Sampling diagnostics (divergences, treedepth, EBFMI)
   - Visual diagnostics (traceplots, densities)
   - Problematic parameters/subjects tables

4. **Recommendations**
   - Actionable next steps
   - Prioritized by severity
   - Model-specific suggestions

5. **Technical Details**
   - Sampling configuration
   - Data filtering parameters
   - Threshold values used
   - Software versions

---

## Understanding Diagnostic Status

### PASS ✓
- All diagnostics within acceptable ranges
- Proceed with confidence
- Results are reliable

### WARN ⚠
- Some diagnostics show warnings
- Review carefully before proceeding
- May be acceptable depending on use case

### FAIL ✗
- Critical diagnostic failures
- DO NOT use for inference
- Address issues before proceeding

---

## Troubleshooting

### Common Issues

**Issue:** "Could not find fit file"
- Check that fit file exists in expected location
- Use `--fit_file` to specify exact path
- Verify task, cohort, model, group names are correct

**Issue:** "Fit type mismatch"
- Ensure `--group` matches actual fit type
- sing = single subject
- batch_XXX = batch of subjects
- hier/group = hierarchical model

**Issue:** "Failed to render HTML"
- Check that rmarkdown is installed
- Check for errors in R console output
- RMD file is still generated even if HTML fails

**Issue:** "Out of memory"
- Large hierarchical models may need more RAM
- Use `--save_rds FALSE` to reduce memory
- Process in smaller batches

### Getting Help

1. Run with `--verbose` for detailed output
2. Check log files in `Outputs/{TASK}/fits/diagnostics/{COHORT}/logs/`
3. Use `--dry_run` to preview what will happen
4. Verify fit object loads correctly in R

---

## Advanced Usage

### Programmatic Access

```r
# Load helpers
source("scripts/diagnostics/helpers/diagnostic_helpers.R")
source("scripts/diagnostics/analysis/single_fit_diagnostics.R")

# Load fit
fit <- readRDS("path/to/fit.rds")

# Run analysis
results <- analyze_single_fit(fit, subject_id = "0001")

# Check status
print(results$overall_status)
print(results$recommendations)

# Access specific diagnostics
print(results$convergence$rhat$max)
print(results$sampling$divergences$rate)
```

### Custom Thresholds

```r
# Create custom thresholds
custom_thresholds <- list(
  thresholds = list(
    rhat = list(problematic = 1.2),  # More lenient
    ess_ratio = list(problematic = 0.2)  # More strict
  )
)

# Use in analysis
results <- analyze_single_fit(fit, thresholds = custom_thresholds)
```

---

## Performance

- **Single fit**: ~10-20 seconds
- **Batch (100 subjects)**: ~30-80 seconds
- **Batch processing (10 models, parallel)**: ~2-5 minutes

---

## Implementation Details

### Architecture

```
Phase 1: Infrastructure (helpers, config)
Phase 2: Analysis (single, batch, hierarchical)
Phase 3: Visualization & Reporting (plots, templates)
Phase 4: Orchestration (run scripts)
```

### Key Design Decisions

- **BIDS naming** - Consistent, parseable filenames
- **Modular design** - Each phase independent
- **Configurable** - YAML-based settings
- **Educational** - Explains diagnostics inline
- **Scalable** - Handles 1 to 1000+ subjects

### Total Implementation

- **24 files created**
- **~6,000 lines of code**
- **4 development phases**
- **Production-ready**

---

## Citation

If you use this pipeline in your research, please cite:

```
MCMC Diagnostics Pipeline for Cognitive Choice Modeling
[Your Project Name]
2025
```

---

## Support

For issues or questions:
1. Check this documentation
2. Review log files
3. Use `--dry_run` and `--verbose`
4. Contact project maintainer

---

## Version History

- **v1.0** (2025-01-09) - Initial complete implementation
  - All four phases implemented
  - Single, batch, and hierarchical support
  - Comprehensive reporting
  - Educational content
  - Batch processing
  - Complete documentation

---

**Pipeline Status: ✓ Complete and Production-Ready**
