# Empirical Bayes Pipeline

This directory contains scripts for implementing empirical Bayes estimation of cognitive model parameters, following the methodology described in Weigard et al. (2021).

## Overview

The pipeline consists of 4 steps:

1. **Subject Selection** - Select a subset of participants for hierarchical prior generation
2. **Hierarchical Fitting** - Fit hierarchical model to the selected subset
3. **Prior Generation** - Extract posterior distributions and calculate informative priors
4. **Empirical Bayes Fitting** - Fit individual models using informative priors

## Scripts

### Individual Step Scripts

- `select_empbayes_subjects.R` - Step 1: Select subjects after data quality filtering
- `fit_empbayes_hierarchical.R` - Step 2: Fit hierarchical model to subset
- `generate_empbayes_priors.R` - Step 3: Generate informative priors from posteriors
- `fit_empbayes_batch.R` - Step 4: Fit individual subjects with priors

### Master Orchestrator

- `run_empbayes_pipeline.R` - Run all steps or individual steps in sequence

## Usage

### Run Complete Pipeline

```bash
Rscript scripts/empbayes/run_empbayes_pipeline.R \
  --step all \
  --task igt_mod \
  --model ev_ddm \
  --source ahrb \
  --n_hier 300 \
  --subjects all \
  --seed 12345
```

### Run Individual Steps

```bash
# Step 1 only
Rscript scripts/empbayes/run_empbayes_pipeline.R \
  --step 1 \
  --task igt_mod \
  --source ahrb \
  --n_hier 300

# Steps 2-4
Rscript scripts/empbayes/run_empbayes_pipeline.R \
  --step 2-4 \
  --task igt_mod \
  --model ev_ddm \
  --source ahrb
```

### Run Scripts Directly

```bash
# Step 1: Select subjects
Rscript scripts/empbayes/select_empbayes_subjects.R \
  -k igt_mod \
  -s ahrb \
  --n_hier 300 \
  --seed 12345

# Step 2: Fit hierarchical model
Rscript scripts/empbayes/fit_empbayes_hierarchical.R \
  -k igt_mod \
  -m ev_ddm \
  -s ahrb

# Step 3: Generate priors
Rscript scripts/empbayes/generate_empbayes_priors.R \
  -k igt_mod \
  -m ev_ddm \
  -s ahrb

# Step 4: Fit with empirical Bayes
Rscript scripts/empbayes/fit_empbayes_batch.R \
  -k igt_mod \
  -m ev_ddm \
  -s ahrb \
  --subjects all \
  --parallel \
  --cores 8
```

## Key Parameters

### Subject Selection
- `--n_hier` - Number of subjects for hierarchical model (e.g., 300)
- `--hier_subs_file` - Use pre-specified subject list instead of random sampling
- `--seed` - Random seed for reproducibility

### Data Quality Filtering
- `--n_trials` - Minimum number of trials required (default: 120)
- `--RTbound_min_ms` - RT minimum bound in milliseconds (default: 50)
- `--RTbound_max_ms` - RT maximum bound in milliseconds (default: 2500)
- `--rt_method` - RT preprocessing method (default: "remove")

### Fitting Parameters
- `--fitting_method` - Fitting method: "mcmc" or "pathfinder" (default: mcmc)
  - **MCMC**: Full Bayesian inference with checkpointing
  - **Pathfinder**: Fast variational approximation using L-BFGS + PSIS

#### MCMC-specific Parameters
- `--n_warmup` - Warmup iterations (default: 3000)
- `--n_iter` - Total iterations (default: 15000)
- `--n_chains` - Number of chains (default: 4)
- `--adapt_delta` - Adapt delta (default: 0.95)
- `--max_treedepth` - Max tree depth (default: 12)

#### Pathfinder-specific Parameters
- `--pf_num_paths` - Number of independent pathfinder runs (default: 4)
- `--pf_draws` - Final draws after PSIS resampling (default: 1000)
- `--pf_single_path_draws` - Draws per single path (default: 250)

### Empirical Bayes Fitting
- `--subjects` - Subjects to fit: "all", "1-100", "1,5,10-20"
- `--n_subs` - Fit first N eligible subjects
- `--parallel` - Use parallel processing
- `--cores` - Number of cores for parallel (default: 4)

## Output Structure

```
Outputs/TASK/empbayes/
├── subs/
│   └── task-{task}_cohort-{cohort}_ses-{ses}_desc-hier_subs.txt
├── hierarchical/
│   └── task-{task}_cohort-{cohort}_ses-{ses}_group-hier_model-{model}_type-fit_desc-output.rds
├── priors/
│   ├── task-{task}_cohort-{cohort}_ses-{ses}_model-{model}_priors.csv
│   └── task-{task}_cohort-{cohort}_ses-{ses}_model-{model}_priors.rds
└── individual/
    └── [temporary individual fits - deleted after combining]

Outputs/TASK/fits/fit/{cohort}/ses-{ses}/
└── task-{task}_cohort-{cohort}_ses-{ses}_group-emp_model-{model}_type-fit_desc-output.rds
```

## Fitting Methods: MCMC vs Pathfinder

### When to Use MCMC
- **Best for:** Final publication-quality estimates
- **Pros:** Full Bayesian inference, accurate uncertainty quantification, convergence diagnostics
- **Cons:** Can be very slow for complex hierarchical models
- **Use when:** You need the most accurate estimates and have sufficient computational time

### When to Use Pathfinder  
- **Best for:** Fast prototyping and development with complex models
- **Pros:** Orders of magnitude faster than MCMC, still captures parameter correlations via importance sampling
- **Cons:** Approximate inference, may underestimate uncertainty slightly
- **Use when:** Hierarchical MCMC is too slow/infeasible, you need quick turnaround for testing
- **Note:** For empirical Bayes, Pathfinder priors are **much better than individual-only estimates**
- **Diagnostics:** Pathfinder provides ELBO (Evidence Lower Bound) values and log probability, stored in the fit object

### Recommendation
For complex models where hierarchical MCMC is impractical:
1. Use Pathfinder for Step 2 (hierarchical prior generation)
2. Priors will still provide substantial shrinkage and information pooling
3. Step 4 (individual fits) still uses MCMC for final estimates

### Pathfinder Output
Pathfinder fits include:
- **`draws`**: Posterior samples (after PSIS resampling)
- **`sampler_diagnostics`**: Contains `lp__` (log probability) and `lp_approx__` (ELBO) draws
- **`diagnostic_summary`**: Summary statistics including:
  - `elbo_mean`, `elbo_sd`: ELBO statistics (quality of approximation)
  - `lp_mean`, `lp_sd`: Log probability statistics
  - `num_paths`: Number of independent pathfinder runs
- **`diagnostics`**: Per-parameter rhat, ESS (same as MCMC)

## Important Notes

### Data Quality Filtering
- The same data quality filters must be applied in Step 1 (subject selection) and Step 4 (fitting)
- This ensures the subjects selected for the hierarchical model are actually eligible for fitting
- Filters are applied BEFORE random sampling to avoid selecting subjects that will be excluded later

### Parameter Scale
- This pipeline assumes models use unconstrained (`_pr`) parameters
- Priors are calculated as simple mean and SD of pooled posteriors
- No truncated normal fitting required for unbounded parameters

### Subject Lists
- Step 1 saves the selected subjects to a text file
- Step 2 loads this file and saves the subject list in the fit object
- This ensures traceability of which subjects were used for prior generation

### Moving Fits
- Individual empirical Bayes fits are initially saved to `empbayes/individual/`
- After combining, the batch file is MOVED (not copied) to the standard `fits/` directory
- Individual files are deleted after successful combination

## Dry Run Mode

All scripts support `--dry_run` to preview what would happen without actually processing data or fitting models:

```bash
Rscript scripts/empbayes/run_empbayes_pipeline.R \
  --step all \
  --task igt_mod \
  --model ev_ddm \
  --source ahrb \
  --n_hier 300 \
  --dry_run
```

## Example Workflow

### Example 1: MCMC (Full Bayesian)

```bash
# Run complete pipeline with 300 subjects for hierarchical, fit all eligible subjects
Rscript scripts/empbayes/run_empbayes_pipeline.R \
  --task igt_mod \
  --model ev_ddm \
  --source ahrb \
  --n_hier 300 \
  --subjects all \
  --seed 12345 \
  --parallel \
  --cores 16
```

### Example 2: Pathfinder (Fast Approximation)

```bash
# Run pipeline using Pathfinder for hierarchical model (much faster)
Rscript scripts/empbayes/run_empbayes_pipeline.R \
  --task igt_mod \
  --model ev_ddm \
  --source ahrb \
  --n_hier 300 \
  --subjects all \
  --seed 12345 \
  --fitting_method pathfinder \
  --pf_num_paths 4 \
  --pf_draws 1000 \
  --parallel \
  --cores 16

# 2. Or run steps separately for more control
Rscript scripts/empbayes/select_empbayes_subjects.R \
  -k igt_mod -s ahrb --n_hier 300 --seed 12345

Rscript scripts/empbayes/fit_empbayes_hierarchical.R \
  -k igt_mod -m ev_ddm -s ahrb

Rscript scripts/empbayes/generate_empbayes_priors.R \
  -k igt_mod -m ev_ddm -s ahrb

Rscript scripts/empbayes/fit_empbayes_batch.R \
  -k igt_mod -m ev_ddm -s ahrb \
  --subjects all --parallel --cores 16
```

## References

Weigard, A., Beltz, A., Reddy, S. N., Wilson, S., & Sripada, C. (2021). Cognitive process modeling addresses context independence violations in the ABCD Study stop-signal task. *Developmental Cognitive Neuroscience*, 50, 100977.
