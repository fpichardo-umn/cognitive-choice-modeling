# Cognitive Choice Modeling Pipeline for Decision-Making Tasks

### Work in Progress

Please note that this pipeline and its documentation are currently under active development. The structure and functionality may change, and some sections of this README may be incomplete. Feedback and contributions are welcome!

### Description
This repository provides a flexible and reproducible framework for the computational modeling of behavioral data from decision-making tasks, such as the Iowa Gambling Task (IGT). ***While originally developed specifically for the modified IGT (hence the repository's original name), it has since been generalized to be more broadly applicable.*** The pipeline is built using R and Stan (via CmdStanR) and includes modules for model fitting, validation through parameter recovery, and posterior predictive checks.

## Features

- **Modular Pipelines:** Separate, script-based workflows for model fitting, parameter recovery, and posterior predictive checks.
- **Bayesian Modeling:** Leverages the power of Stan for robust Bayesian analysis.
- **Model Validation:** Includes built-in workflows for parameter recovery and posterior predictive checks to ensure model validity.
- **HPC-Ready:** Integrated with SLURM submission scripts for efficiently running large-scale analyses on High-Performance Computing (HPC) clusters.
- **Organized Structure:** Aims for a BIDS-inspired data and analysis structure to enhance reproducibility and clarity (work-in-progress).

## Repository Structure

The repository is organized to separate model code, scripts, and data outputs. Here are the most important directories:

```
.
├── analysis/         # Output location for rendered reports (.Rmd, .html)
├── data/             # Raw task data not included in this repository but expected
│   └── raw/
├── outputs/             # Default output location for generated data from pipelines
│   ├── igt/
│   └── igt_mod/             
├── models/TASK/STATUS      # Stan model definitions organized by task and status of the model: canonical, experimental, working (not tracked)
│   ├── bin/                # Compiled Stan models
│   └── txt/                # Human-readable .stan model files
└── scripts/                # All operational scripts (R and shell)
│   ├── compile/            # Scripts to compile Stan models
│   ├── configs/            # Configuration files for runs (e.g., MCMC settings)
│   ├── fit/                # Scripts for fitting models to data
│   ├── helpers/            # Utility functions used across pipelines
│   ├── model_comparison/   # Scripts to aggregate and compare results
│   ├── parameter_recovery/ # Scripts for the parameter recovery workflow
│   ├── ppc/                # Scripts for posterior predictive checks
└── └── submit_scripts/     # SLURM job submission scripts
```

## Prerequisites & Setup

To use this pipeline, you will need R, CmdStanR, and the R packages listed in the `DESCRIPTION` file.

1.  **Clone the Repository:**
    ```bash
    git clone [https://github.com/fpichardo-umn/cognitive-choice-modeling.git](https://github.com/fpichardo-umn/cognitive-choice-modeling.git)
    cd cognitive-choice-modeling
    ```

2.  **Install R and CmdStanR:**
    Follow the official instructions to get [R](https://www.r-project.org/) and [CmdStanR](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) installed and configured on your system.

3.  **Install R Packages:**
    Open an R session in the project's root directory and install the required packages.
    ```r
    # In your R console
    install.packages(c("optparse", "here", "data.table", "ggplot2", "patchwork", "ggridges", "bayesplot", "posterior", "loo", "stringr", "fs", "purrr", "tidyr", "dplyr", "readr"))
    ```

4.  **Compile the Stan Models:**
    The Stan models must be compiled before they can be used. Run the following script from your terminal, specifying the task name.
    ```bash
    Rscript scripts/compile/compile_models_all.R --task igt_mod
    ```

## Data Preparation

The pipeline expects raw data to be placed in the `data/raw/{cohort_name}/` directory, though this is not tracked by Git. The processing scripts expect data files in a `.csv` format with the following essential columns:

-   `sub_id`: A unique subject identifier.
-   `trial`: The trial number, starting from 1.
-   `choice`: The deck or option chosen (e.g., 1, 2, 3, 4).
-   `outcome`: The reward or loss received on that trial.

Additional columns for covariates can be included and specified in the modeling scripts.

## Core Workflows (Usage)

The pipeline is designed to be run from the command line, making it ideal for scripting and HPC use. Below are examples for the primary workflows.

* **Fitting Workflow:** Takes raw data and produces fitted model parameters.
* **Validation Workflows:** Use fitted models (for Parameter Recovery and PPC) or simulated data (for Parameter Recovery) to ensure the models are behaving as expected.

#### Workflow 1: Fitting a Model

To fit a single model to a single subject's data, use the `fit_single_model_single_subject.R` script.

**Example:**
```bash
Rscript scripts/fit/fit_single_model_single_subject.R \
  --model "vpp_reg" \
  --task "igt_mod" \
  --sub "1001" \
  --ses "01" \
  --fit_path "data"
```
For large-scale fitting on an HPC, use the corresponding SLURM submission script:

```Bash
sbatch scripts/submit_scripts/sbatch_fit.sh
```
(You will need to edit the script to configure the array of subjects, models, etc.)

#### Workflow 2: Parameter Recovery
Parameter recovery is a multi-step process to validate your models. A wrapper script, `run_full_PR.sh`, helps automate this. The general process is:

1. Generate known parameters (`01_generate_parameters.R`)
2. Simulate behavioral data from those parameters (`02_simulate_data.R`).
3. Fit the models to the simulated data (`03_run_fitting.sh`).
4. Summarize and compare the recovered vs. true parameters (`04_summarize_recovery.R`).

To run the full pipeline on an HPC, use the SLURM submission script:

```Bash
sbatch scripts/submit_scripts/sbatch_pr.sh
```
#### Workflow 3: Posterior Predictive Checks (PPC)
PPC helps validate a model by checking if it can generate data that looks like the real data it was fit to. This also calculates the log-likelihood for model comparison (e.g., LOOIC).

**Example:**

```Bash
Rscript scripts/ppc/ppc_and_ll_single_model_single_subject.R \
  --model "vpp_reg" \
  --task "igt_mod" \
  --sub "1001" \
  --ses "01" \
  --fit_path "data" \
  --ppc_path "data"
```
For large-scale PPC on an HPC, use the SLURM submission script:

```Bash
sbatch scripts/submit_scripts/sbatch_ppc.sh
```
#### Model Comparison
After running the PPC pipeline for multiple models, you can aggregate the `loo` results to compare model fit using the `gather_and_compare_models_loo.R` script.

**Example:**

```Bash
Rscript scripts/model_comparison/gather_and_compare_models_loo.R --task "igt_mod"
```
This will generate a summary CSV file comparing the models in the `analysis/model_comparison` directory.
