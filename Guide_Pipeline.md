# Computational Modeling Pipeline User Guide

*For researchers familiar with computational modeling who want to learn this pipeline*

---

## Pipeline Overview

This pipeline provides a complete research infrastructure for cognitive modeling with three core capabilities:

1. **Model Fitting:** Estimate parameters for several cognitive models - or you can build your own
2. **Model Validation:** Parameter recovery and posterior predictive checks  
3. **Model Comparison:** Systematic comparison across multiple criteria

### IMPORTANT NOTE:
- This guide was generated primarily via a Claude chat
- Please review the required arguments for each script prior to running

**Key Features:**
- **Scale:** Individual subjects → HPC clusters
- **Validation:** Built-in parameter recovery and behavioral checks
- **Organization:** BIDS-compliant file structure and naming
- **Reproducibility:** Configuration management and automated workflows

**Important Note:** Integrating new datasets and models requires specific data formatting and may need pipeline modifications. The system is designed for specific experimental paradigms and data structures.

---

## Quick Start

### Prerequisites
- **Data format:** Your data should include `subjID`, `trial`, `choice`, `wins`, `losses` columns
- **Environment:** R with cmdstanr, Stan models compiled for your system
- **Computing:** Can run locally or on SLURM clusters

### Your First Model Fit
```bash
# Fit Expected Value model to one subject
Rscript scripts/fit/fit_single_model.R \
  -m ev \                          # Model: Expected Value
  -k igt \                         # Task: Iowa Gambling Task  
  -s your_cohort \                 # Your data identifier
  --subid 1001                     # Subject ID

# Results: Data/igt/rds/fit/your_cohort/task-igt_cohort-your_cohort_group-sing_model-ev_sub-1001_type-fit.rds
```

### Validate the Model
```bash
# Test parameter recovery before trusting results
./scripts/parameter_recovery/run_full_PR.sh \
  -k igt \                             # Task
  -m ev \                              # Model to test
  -s your_cohort \                     # Data source for realism
  -g sing \                            # Individual fitting approach
  pr_genparams,pr_simulate,pr_recovery # All steps: generate → simulate & fit → analyze

# Results: Data/igt/sim/recovery/ with correlation/bias statistics
```

### Check Behavioral Realism
```bash
# Posterior predictive checks
./scripts/ppc/run_ppc.sh \
  --model ev \                     # Your fitted model
  --cohort your_cohort \           # Data source
  --task igt \                     # Task type
  --steps all                      # Full PPC pipeline

# Results: Data/ppc/reports/ with behavioral validation
```

---

## Understanding the Model Library

### Directory Structure
```
models/
├── igt/                        # Iowa Gambling Task models
│   ├── canonical               # Status of model: canonical, experimental, and working (not tracked)
│     ├── txt/                  # Stan source files (.stan)
│     │   └── fit/              # Main model files
│     └── bin/                  # Compiled executables (system-specific)
│         └── fit/
└── igt_mod/                    # Modified IGT models (play/pass paradigm)
```

### Model Naming Convention
```
task-igt_group-sing_model-pvldelta_type-fit.stan
│    │   │     │           │         │        │
│    │   │     │           │         │        └─ File extension
│    │   │     │           │         └────── Model type (fit/postpc/prepc)
│    │   │     │           └─────────────── Model name
│    │   │     └─────────────────────── Fitting approach (sing/batch/hier)
│    │   └──────────────────────────── Group type
│    └────────────────────────────────── Task name
└────────────────────────────────────── BIDS prefix
```

### Available Models
**Model Statuses:**
It's beneficial to organize the models by their validation status to maintain a better set of models and work on others without uploading to the official pipeline. If you know the status, you can pass it where needed. However, the pipeline is able to auto-detect which status the specific model is in if it's been compiled.
- **Canonical**: These are models that have been used previously in the literature
- **Experimental**: These are models that are being tested with this pipeline but have not been published
- **Working**: This is a working directory of models that are not tracked or maintained on the official pipeline but may be referenced in some code

**Reinforcement Learning Canonical Models:**
- **`ev`** - Expected Value: Basic utility tracking
- **`pvldelta`** - Prospect Valuation Learning with delta rule
- **`pvldecay`** - PVL with decay learning  
- **`pvlboth`** - PVL with delta + decay
- **`orl`** - Outcome Representation Learning: Separate reward/punishment learning systems
- **`vpp`** - Value-Plus-Perseveration: Adds choice perseveration
- **`vse`** - Value and Sequential Exploration: Explicit exploration

**Sequential Sampling Canonical Models:**
- **`ddm`** - Drift Diffusion Model: Models reaction times

**Hybrid Experimental Models (RL + SSM):**
- **`ev_ddm`**, **`pvldelta_ddm`**, **`orl_ddm`**, etc.
- Currently developed for `igt_mod`, under development for `igt`

**Fitting Approaches:**
- **`sing`** - Individual models (one per subject)
- **`hier`** - Hierarchical models (population + individual parameters)

### Finding Models
```bash
# List all available models for IGT
ls models/igt/*/txt/fit/

# Find models with specific features
ls models/igt/*/txt/fit/ | grep "pvl"     # All PVL variants
ls models/igt/*/txt/fit/ | grep "hier"    # All hierarchical variants

# Check if model is compiled
ls models/igt/*/bin/fit/task-igt_group-sing_model-ev_type-fit
```

---

## Core Workflows

### Workflow 1: Model Fitting

#### Individual Subject Fitting
```bash
# Detailed analysis of single subjects
Rscript scripts/fit/fit_single_model.R \
  -m pvldelta \                    # Model name
  -k igt \                         # Task
  -s cohort_name \                 # Data source
  --ses session \                  # Session number (00, 01)
  --subid 1001 \                   # Subject ID
  --model_status canonical \       # Model validation status - this is not required!
  --n_warmup 1000 \                # MCMC warmup
  --n_iter 2000 \                  # MCMC iterations
  --n_chains 4                     # MCMC chains

# Advanced options:
  --adapt_delta 0.95 \             # Step size adaptation
  --max_treedepth 12 \             # Maximum tree depth
  --seed 12345                     # Reproducibility
```

#### Batch Fitting (Many Subjects)
```bash
# Submit parallel jobs for multiple subjects
./scripts/submit_scripts/submit_batch_models.sh \
  -a "1-100" \                     # Subject range or list
  -m pvldelta \                    # Model name
  -k igt \                         # Task
  -s cohort_name \                 # Data source
  -S session \                     # Session number
  -e your@email.edu                # Email for notifications

# Advanced batch options:
  -f default \                     # Configuration profile
  -n                               # Dry run (show commands without executing)
```

#### Hierarchical Fitting
```bash
# Population + individual parameters simultaneously
Rscript scripts/fit/fit_hierarchical_models_cmdSR.R \
  -m pvldelta \                    # Model
  -k igt \                         # Task
  -s cohort_name \                 # Data source
  --ses session \                  # Session number (00, 01)
  --n_subs 100 \                   # Number of subjects
  --n_trials 100                   # Trials per subject

# Hierarchical-specific options:
  --group hier \                   # Hierarchical group type
  --check_iter 1000                # Checkpoint interval
```

### Workflow 2: Parameter Recovery Validation

**Complete 5-step process:** Fit real data → Generate parameters → Simulate data → Fit simulated data → Analyze recovery

```bash
# Full parameter recovery pipeline (starts with fitting your real data)
./scripts/parameter_recovery/run_full_PR.sh \
  -k igt \                         # Task
  -m pvldelta \                    # Model to test
  -s cohort_name \                 # Data source
  --ses session \                  # Session number (00, 01)
  -g sing \                        # Fitting approach (sing/hier)
  all                              # Run all steps: fit → pr_gen → pr_sim → pr_fit → pr_recovery

# Individual steps (advanced usage):
./scripts/parameter_recovery/run_full_PR.sh -k igt -m pvldelta -s cohort --ses session -g sing fit         # Fit real data
./scripts/parameter_recovery/run_full_PR.sh -k igt -m pvldelta -s cohort --ses session -g sing pr_gen      # Generate parameters
./scripts/parameter_recovery/run_full_PR.sh -k igt -m pvldelta -s cohort --ses session -g sing pr_sim      # Simulate data
./scripts/parameter_recovery/run_full_PR.sh -k igt -m pvldelta -s cohort --ses session -g sing pr_fit      # Fit simulated data
./scripts/parameter_recovery/run_full_PR.sh -k igt -m pvldelta -s cohort --ses session -g sing pr_recovery # Analyze recovery

# Just fit real data (no parameter recovery validation):
./scripts/parameter_recovery/run_full_PR.sh -k igt -m pvldelta -s cohort -g sing fit
```

**Interpreting recovery results:**
```r
# Good recovery indicators:
correlation > 0.8     # Strong true vs. recovered correlation
bias < 0.1           # Small systematic error  
rmse < 0.2           # Small root mean square error

# Poor recovery warnings:
correlation < 0.6     # Weak relationship suggests identification problems
bias > 0.3           # Large systematic error suggests model bias
rmse > 0.5           # Large overall error suggests poor precision
```

**Important:** Models may recover some parameters better than others. A model might excel at recovering learning parameters but struggle with motivational parameters. Consider parameter-specific recovery quality when interpreting results.

### Workflow 3: Posterior Predictive Checks

**Test long-term behavioral realism of fitted models:**

```bash
# Complete PPC pipeline
./scripts/ppc/run_ppc.sh \
  --model pvldelta \               # Fitted model
  --cohort cohort_name \           # Data source
  --ses session \                  # Session number (00, 01)
  --task igt \                     # Task type
  --steps all \                    # simulation → stats → loglik → report
  --n_sims 100 \                   # Simulations per subject
  --sampling weighted              # Posterior sampling method

# PPC-specific options:
  --block_size 20 \                # Trials per block for analysis
  --ic_method loo \                # Information criterion (loo/waic)
  --render                         # Generate HTML report
```

**Behavioral statistics calculated:**
```r
# IGT-specific metrics:
deck1_freq, deck2_freq, deck3_freq, deck4_freq    # Deck selection frequencies
good_deck_freq, bad_deck_freq                     # Good vs. bad deck preferences  
win_stay, lose_shift                              # Strategy measures
net_score                                         # Performance measure
perseveration                                     # Choice persistence

# RT metrics (for DDM models):
rt_mean, rt_sd, rt_q10, rt_q50, rt_q90           # Reaction time distributions
```

**Interpreting PPC results:**
```r
# Good behavioral fit:
ppp_values ≈ 0.5       # Posterior predictive p-values around 0.5
extreme_failures < 5%   # Few statistics with ppp < 0.05 or > 0.95

# Poor behavioral fit indicators:
ppp_values far from 0.5 # Systematic over/under-prediction
extreme_failures > 20%  # Many poorly captured statistics
```

### Workflow 4: Model Comparison

**Systematic comparison across multiple criteria:**

```bash
# Compare specific models
./scripts/model_comparison/run_model_comparison.R \
  --task igt \
  --cohort cohort_name \
  --ses session \
  --models "ev,pvldelta,orl" \     # Models to compare
  --analysis_types recovery,ppc,ic,report

# Comprehensive comparison
./scripts/model_comparison/run_model_comparison.R \
  --task igt \
  --cohort cohort_name \
  --ses session \
  --models all \                   # All available models
  --analysis_types recovery,ppc,ic,report \
  --render_html                    # Generate HTML report

# Model comparison options:
  --model_types rl_only,hybrid \   # Filter by model type
  --comparison_name study1 \       # Custom identifier
  --save_plots                     # Save individual plots
```

**Comparison criteria:**
1. **LOOIC/WAIC:** Short-term prediction (predict next choice for this dataset)
2. **Parameter Recovery:** Model validity in general (can it reliably estimate parameters?)
3. **Posterior Predictive Checks:** Long-term prediction (full behavioral patterns for this dataset)

**Key insight:** No single model may be "best" across all criteria. A model might excel at parameter recovery but struggle with behavioral prediction, or vice versa.

---

## Configuration Management

### Configuration Profiles

**Fitting configurations (`scripts/configs/`):**
```bash
fit_params_quick.conf:    n_warmup=100,  n_iter=1000,  n_chains=2    # Fast testing
fit_params_default.conf:  n_warmup=1000, n_iter=2000,  n_chains=4    # Standard
fit_params_sing.conf:     n_warmup=500,  n_iter=2000,  n_chains=4    # Individual models
```

**Data configurations:**
```bash
data_params_default.conf: n_trials=100, RTbound_ms=50, rt_method="remove"
data_params_quick.conf:   n_trials=60,  RTbound_ms=100
```

**Using configurations:**
```bash
# Specify configuration in batch jobs
./submit_batch_models.sh -a "1-10" -m ev -k igt -s cohort --ses session -f quick -d default

# Or in individual scripts
Rscript fit_single_model.R -m ev -k igt --subid 1001 --config quick
```

### File Organization (BIDS-Compliant)

**Input data:**
```
Data/raw/COHORT_NAME/
├── TASK_COHORT_SES.csv          # Your behavioral data - e.g., igt_es_04.csv
└── ses-SES/subject_ids_all.txt          # Subject identifier lists
```

**Model fitting results:**
```
Outputs/igt/fits/fit/cohort_name/
├── task-igt_cohort-name_ses-session_group-sing_model-ev_sub-1001_type-fit.rds
├── task-igt_cohort-name_ses-session_group-batch_model-ev_type-fit.rds
└── task-igt_cohort-name_ses-session_group-hier_model-ev_type-fit.rds
```

**Simulation results:**
```
Outputs/igt/simulation/
├── data/rds                      # Simulated behavioral data
├── data/txt                      # Simulated behavioral data
└── parameters/                   # Generated parameter setsRecovery analysis results
```

**Parameter recovery results:**
```
Outputs/igt/validation/parameter_recovery/
├── fits/                       # Estimated fits for the simulated data
└── analysis/                   # Recovery analysis results
```

**PPC results:**
```
Data/ppc/cohort-name/ses-SES
├── sim/                        # Simulated data from fitted models
├── stats/                      # Behavioral statistics
├── loglik/                     # Log-likelihood values  
└── reports/                    # Generated reports
```

---

## Advanced Features

### High-Performance Computing

**SLURM integration:**
```bash
# Automatic job submission with resource management
./submit_batch_models.sh -a "1-500" -m pvldelta -k igt -s cohort --ses session -e you@email.edu

# Resources automatically allocated:
#SBATCH --nodes=1 --ntasks-per-node=8 --mem=32gb --time=24:00:00
#SBATCH --mail-type=ALL --mail-user=you@email.edu
```

**Checkpointing and recovery:**
- Models automatically save progress at specified intervals
- Jobs can resume from last checkpoint if interrupted
- Email notifications for job completion/failure

### Model Compilation

**Intelligent compilation system:**
```bash
# Compile all models for your task
./scripts/compile/compile_models_cmdSR.R -k igt -t all --yes

# Compile specific models with pattern matching
./scripts/compile/compile_models_cmdSR.R -k igt -m "pvl" -s "ddm" -e "test" --yes

# Hash-based change detection (only recompiles when Stan code changes)
# BIDS filename conversion (legacy → modern naming)
```

### Custom Analysis Development

**Building on pipeline infrastructure:**
```r
# Load helper functions
source("scripts/helpers/helper_functions_cmdSR.R")

# Use existing data loading
data <- load_data(task = "igt", cohort = "your_cohort")

# Use existing model registry
model_info <- get_model_defaults("igt")

# Use existing file organization
output_file <- generate_bids_filename(
  task = "igt", model = "custom", cohort = "your_cohort",
  additional_tags = list(analysis = "custom"), ext = "rds"
)
```

---

## Best Practices

### Model Development Workflow
1. **Start simple, build complexity systematically:**
   - **Model complexity:** EV → VPP → ORL → Hybrid models
   - **Fitting approach:** Individual (sing) → Hierarchical (hier)  
   - **Model type:** RL → RL+DDM hybrid
2. **Validate early with Parameter Recovery:** Use PR to eliminate models from consideration since PPC takes much longer
3. **Expect multiple adequate models:** You'll likely have several models that perform well
4. **Understand trade-offs:** Models may excel at recovering some parameters/behaviors better than others - no single "winner"
5. **Use all three criteria:** LOOIC (short-term prediction), PR (general validity), PPC (long-term prediction)

### Troubleshooting Common Issues

**Fitting problems:**
```r
# Check convergence diagnostics
fit$diagnostics$rhat_max        # Should be < 1.1
fit$diagnostics$ess_min         # Should be > 400

# Common solutions:
# - Increase adapt_delta (0.8 → 0.95 → 0.99)
# - Increase max_treedepth (10 → 12 → 15)  
# - Increase n_warmup and n_iter
# - Check for data issues (missing values, outliers)
```

**Parameter recovery failures:**
```r
# Check parameter identifiability
# - Some parameters may be confounded
# - Model may be too complex for data
# - Consider simpler model or more data

# Check simulation accuracy
# - Verify task implementation matches real experiment
# - Check parameter ranges are realistic
```

**PPC failures:**
```r
# Model may be missing key psychological processes
# - Try more complex model family
# - Check if hybrid RL+DDM model needed
# - Consider task-specific modifications
```

### Performance Optimization

**For large datasets:**
- **Hierarchical models:** May be too slow for complex models with large datasets
  - Try hierarchical first to test feasibility
  - Fall back to batch individual fits if hierarchical is too slow
  - Consider two-step or empirical Bayesian hierarchical estimation
  - Consider alternatives to MCMC for very large datasets
- **Batch processing:** Use parallel individual fits when hierarchical is impractical
- **HPC resources:** Essential for large-scale analyses

**For development:**
- Use quick configuration profiles for testing
- Test on subset of subjects before full analysis
- Use dry-run options to check commands before execution

---

## Getting Help

### Documentation
- **Model details:** Check `models/igt/txt/fit/` for Stan code
- **Helper functions:** See `scripts/helpers/` for implementation details
- **Examples:** Look at submit scripts for common usage patterns

### File locations for troubleshooting
- **Log files:** `log_files/` directory with detailed execution logs
- **Configuration:** `scripts/configs/` for parameter settings
- **Results:** Follow BIDS naming convention to locate outputs

### Common commands summary
```bash
# Basic fitting
Rscript scripts/fit/fit_single_model.R -m [MODEL] -k [TASK] -s [COHORT] --ses [SESSION] --subid [ID]

# Validation  
./scripts/parameter_recovery/run_full_PR.sh -k [TASK] -m [MODEL] -s [COHORT] --ses [SESSION] -g sing all

# Behavioral checks
./scripts/ppc/run_ppc.sh --model [MODEL] --cohort [COHORT] --ses [SESSION] --task [TASK] --steps all

# Model comparison
./scripts/model_comparison/run_model_comparison.R --task [TASK] --cohort [COHORT] --ses [SESSION] --models [LIST]
```

---

*This guide covers the practical usage of the computational modeling pipeline. For conceptual background on computational modeling, see the Newcomers Guide.*
