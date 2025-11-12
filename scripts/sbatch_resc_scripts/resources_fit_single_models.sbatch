#!/bin/bash -l
#SBATCH -J ${JOB_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --tmp=32gb
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH -o log_files/${TASK}/fit/%x_%A_%4a.out
#SBATCH -e log_files/${TASK}/fit/%x_%A_%4a.err

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

# Get the subject ID based on the array task ID
SUBJECT_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SUBS_LIST_FILE})

if [ -z "$SUBJECT_ID" ]; then
  echo "Error: Unable to get subject ID for array task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

# Determine RT handling
if [ -z "${rt_method}" ]; then
  if [[ "${MODEL_NAME}" == *"ddm"* ]]; then
    rt_method="remove"
  else
    rt_method="force"
  fi
fi

# Run the R script with the configuration
Rscript ./scripts/fit/fit_single_model.R \
  -m ${MODEL_NAME} \
  -t ${MODEL_TYPE} \
  -k ${TASK} \
  -s ${SOURCE} \
  $([ ! -z "${SESSION}" ] && echo "--ses ${SESSION}") \
  --subid ${SUBJECT_ID} \
  --index ${SLURM_ARRAY_TASK_ID} \
  --n_trials ${n_trials} \
  --RTbound_min_ms ${RTbound_min_ms} \
  --RTbound_max_ms ${RTbound_max_ms} \
  --rt_method ${rt_method} \
  --n_warmup ${n_warmup} \
  --n_iter ${n_iter} \
  --n_chains ${n_chains} \
  --adapt_delta ${adapt_delta} \
  --max_treedepth ${max_treedepth} \
  --check_iter ${CHECK_ITER} \
  --seed ${SLURM_JOB_ID} \
  $([ "${enable_adaptive_iter}" = "TRUE" ] && [ ! -z "${min_iter}" ] && echo "--min_iter ${min_iter}") \
  $([ "${enable_adaptive_iter}" = "TRUE" ] && [ ! -z "${max_iter}" ] && echo "--max_iter ${max_iter}") \
  $([ "${enable_adaptive_iter}" = "TRUE" ] && [ ! -z "${iter_increment}" ] && echo "--iter_increment ${iter_increment}") \
  $([ "${enable_adaptive_iter}" = "TRUE" ] && [ ! -z "${target_rhat}" ] && echo "--target_rhat ${target_rhat}") \
  $([ "${enable_adaptive_iter}" = "TRUE" ] && [ ! -z "${target_ess_bulk}" ] && echo "--target_ess_bulk ${target_ess_bulk}") \
  $([ "${enable_adaptive_iter}" = "TRUE" ] && [ ! -z "${target_ess_tail}" ] && echo "--target_ess_tail ${target_ess_tail}") \
  $([ "${enable_adaptive_iter}" = "FALSE" ] && echo "--disable_adaptive_iter") \
  --init