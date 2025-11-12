#!/bin/bash -l
#SBATCH -J ${JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --tmp=32gb
#SBATCH -t 96:00:00
#SBATCH --mail-type=ALL
#SBATCH -o log_files/fit/%x_%A_%a.out
#SBATCH -e log_files/fit/%x_%A_%a.err

# Check if email is provided
if [ -z "$USER_EMAIL" ]; then
  echo "Error: USER_EMAIL is not set. Please provide an email address."
  exit 1
fi

#SBATCH --mail-user=$USER_EMAIL

# Load necessary modules
module load R/4.2.0-rocky8
module load gcc/13.1.0-mptekim
module load cmake/3.29.2-rocky8 
module load boost/1.82.0-gcc-13.1.0-sox6zlc

# Change to the submit directory
cd ${SLURM_SUBMIT_DIR}

# Source the config files
source ./scripts/configs/fit_params_${FIT_CONFIG}.conf
source ./scripts/configs/data_params_${DATA_CONFIG}.conf

# Determine RT handling
if [ -z "${RT_METHOD}" ]; then
  if [[ "${MODEL_NAME}" == *"ddm"* ]]; then
    RT_METHOD="mark"
  else
    RT_METHOD="force"
  fi

fi

# Run the R script with the configuration
Rscript ./scripts/fit/fit_hierarchical_models_cmdSR.R \
  -m ${MODEL_NAME} \
  -t ${MODEL_TYPE} \
  -k ${TASK} \
  -g ${GROUP_TYPE} \
  --n_subs ${n_subs} \
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
    $([ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ] && [ ! -z "${MIN_ITER}" ] && echo "--min_iter ${MIN_ITER}") \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ] && [ ! -z "${MAX_ITER}" ] && echo "--max_iter ${MAX_ITER}") \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ] && [ ! -z "${ITER_INCREMENT}" ] && echo "--iter_increment ${ITER_INCREMENT}") \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ] && [ ! -z "${TARGET_RHAT}" ] && echo "--target_rhat ${TARGET_RHAT}") \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ] && [ ! -z "${TARGET_ESS_BULK}" ] && echo "--target_ess_bulk ${TARGET_ESS_BULK}") \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "TRUE" ] && [ ! -z "${TARGET_ESS_TAIL}" ] && echo "--target_ess_tail ${TARGET_ESS_TAIL}") \
    $([ "${ENABLE_ADAPTIVE_ITER}" = "FALSE" ] && echo "--disable_adaptive_iter")
