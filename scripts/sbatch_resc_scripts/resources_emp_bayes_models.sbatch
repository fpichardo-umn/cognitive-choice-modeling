#!/bin/bash -l
#SBATCH -J ${JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --tmp=32gb
#SBATCH -t 96:00:00
#SBATCH --mail-type=ALL
#SBATCH -o log_files/empbayes/%x_%A_%a.out
#SBATCH -e log_files/empbayes/%x_%A_%a.err

# Check if email is provided
if [ -z "$USER_EMAIL" ]; then
  echo "Error: USER_EMAIL is not set. Please provide an email address."
  exit 1
fi

#SBATCH --mail-user=$USER_EMAIL

# Convert STEPS string to array
IFS=',' read -ra STEP_ARRAY <<< "$STEPS"

# Function to check if a step should be run
should_run_step() {
  for step in "${STEP_ARRAY[@]}"; do
    if [ "$step" = "$1" ]; then
      return 0  # Return true (0 in bash) if the step is found
    fi
  done
  return 1  # Return false (1 in bash) if the step is not found
}

# Load necessary modules
module load R/4.2.0-rocky8
module load gcc/13.1.0-mptekim
module load cmake/3.29.2-rocky8 
module load boost/1.82.0-gcc-13.1.0-sox6zlc

# Change to the submit directory
cd ${SLURM_SUBMIT_DIR}

# Source the config files
source ./scripts/configs/fit_params_emp_hier_${FIT_CONFIG}.conf
source ./scripts/configs/data_params_${DATA_CONFIG}.conf

# Function to check exit status and log
check_status() {
    if [ $1 -ne 0 ]; then
        echo "Error: $2 failed with exit code $1"
        exit $1
    else
        echo "$2 completed successfully"
    fi
}

# Step 1: Run initial hierarchical model
if should_run_step 1; then
  echo "Starting initial hierarchical model fit..."
  Rscript ./scripts/empbayes/fit_empbayes_hierarchical.R \
    -m ${MODEL_NAME} \
    -s ${SOURCE} \
    --ses ${SES} \
    -k ${TASK} \
    --n_trials ${N_TRIALS} \
    --RTbound_min_ms ${RTBOUND_MIN_MS} \
    --RTbound_max_ms ${RTBOUND_MAX_MS} \
    --rt_method ${RT_METHOD} \
    --n_warmup ${N_WARMUP} \
    --n_iter ${N_ITER} \
    --n_chains ${N_CHAINS} \
    --adapt_delta ${ADAPT_DELTA} \
    --max_treedepth ${MAX_TREEDEPTH} \
    --check_iter ${CHECK_ITER} \
    --seed ${SLURM_JOB_ID}
  check_status $? "Initial hierarchical model fit"
fi

# Step 2: Generate informative priors
if should_run_step 2; then
  echo "Generating informative priors..."
  Rscript ./scripts/empbayes/generate_empbayes_priors.R \
    -m ${MODEL_NAME} \
    -s ${SOURCE} \
    --ses ${SES} \
    -k ${TASK}
  check_status $? "Generating informative priors"
fi

# Step 3: Run full empirical Bayes model
if should_run_step 3; then
  echo "Starting full empirical Bayes model fit..."
  # Source the config files
  source ./scripts/configs/fit_params_emp_indiv_${FIT_CONFIG}.conf
  
  # max number of subjects: a
  Rscript ./scripts/empbayes/fit_empbayes_batch.R \
    -m ${MODEL_NAME} \
    -s ${SOURCE} \
    --ses ${SES} \
    -k ${TASK} \
    -a "all" \
    --n_trials ${N_TRIALS} \
    --RTbound_min_ms ${RTBOUND_MIN_MS} \
    --RTbound_max_ms ${RTBOUND_MAX_MS} \
    --rt_method ${RT_METHOD} \
    --n_warmup ${N_WARMUP} \
    --n_iter ${N_ITER} \
    --n_chains ${N_CHAINS} \
    --adapt_delta ${ADAPT_DELTA} \
    --max_treedepth ${MAX_TREEDEPTH} \
    --check_iter ${CHECK_ITER} \
    --subs_file ${SUBS_FILE} \
    --seed ${SLURM_JOB_ID}
  check_status $? "Full empirical Bayes model fit"
fi
