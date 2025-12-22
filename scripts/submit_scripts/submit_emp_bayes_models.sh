#!/bin/bash

# Load modules
module load R/4.2.0-rocky8

# Function to print usage
print_usage() {
  echo "Usage: $0 -m <model_names> -f <fit_config> -d <data_config> -e <email> -k <task> -[-n]"
  echo "Example: $0 -m \"ev,pvl\" -f default -d full -e your@email.edu -k igt_mod"
  echo "Options:"
  echo "  -m    Comma-separated list of model names"
  echo "  -s    Data source (cohort)"
  echo "  -S    Session identifier"
  echo "  -f    Fit parameters config name for empirical bayes (default: simple) [complex]"
  echo "  -d    Data parameters config name (default: default)"
  echo "  -k    Task name (e.g., igt_mod)"
  echo "  -c    Number of iterations for checkpoint runs (default: 10000)"
  echo "  -l    Subs list file [Data/raw/COHORT/ses-SES/] (default: subject_ids_complete_valid.txt)"
  echo "  -t    Steps to run (default: "1,2,3", options: "1", "2", "3", "1,2", "2,3", or any combination)"
  echo "  -e    Your email address (required)"
  echo "  -n    Dry run (optional)"
  exit 1
}

# Parse command line arguments
DRY_RUN=false
while getopts ":m:s:S:f:d:k:c:l:t:e:n" opt; do
  case $opt in
    m) MODEL_NAMES=$OPTARG ;;
    s) SOURCE=$OPTARG ;;
    S) SES=$OPTARG ;;
    f) FIT_CONFIG=$OPTARG ;;
    d) DATA_CONFIG=$OPTARG ;;
    k) TASK=$OPTARG ;;
    c) CHECK_ITER=$OPTARG ;;
    l) SUBS_FILE=$OPTARG;;
    t) STEPS=$OPTARG ;;
    e) USER_EMAIL=$OPTARG ;;
    n) DRY_RUN=true ;;
    \?) echo "Invalid option -$OPTARG" >&2; print_usage ;;
  esac
done

# Check if required arguments are provided
GROUP_TYPE="group_hier" # It's always hier and then emp (which is generated in the scripts)
if [ -z "$MODEL_NAMES" ] || [ -z "SOURCE" ] ||  [ -z "SES" ] ||  [ -z "$USER_EMAIL" ] || [ -z "$TASK" ]; then
  echo "Error: Model names, cohort/source, session nubmer, email address, and task are required."
  print_usage
fi

# Set default values if not provided
FIT_CONFIG=${FIT_CONFIG:-simple}
DATA_CONFIG=${DATA_CONFIG:-default}
MODEL_TYPE=${MODEL_TYPE:-fit}
CHECK_ITER=${CHECK_ITER:-10000}
STEPS=${STEPS:-"1,2,3"}
SUBS_FILE=${SUBS_FILE:-"subject_ids_complete_valid.txt"}

if ! [[ $STEPS =~ ^[1-3](,[1-3]){0,2}$ ]]; then
  echo "Error: Invalid steps specified. Use a comma-separated list of 1, 2, and/or 3."
  exit 1
fi

# Directory checks
SUBMIT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
SCRIPT_DIR="${SUBMIT_DIR}/.."
CONFIG_DIR="${SCRIPT_DIR}/configs"
R_SCRIPT="${SCRIPT_DIR}/empbayes/fit_empbayes_hierarchical.R"
R_SCRIPT2="${SCRIPT_DIR}/empbayes/generate_empbayes_priors.R"
R_SCRIPT3="${SCRIPT_DIR}/empbayes/fit_empbayes_batch.R"
SBATCH_SCRIPTS_DIR="${SCRIPT_DIR}/sbatch_resc_scripts"
SBATCH_SCRIPT="${SBATCH_SCRIPTS_DIR}/resources_emp_bayes_models.sbatch"

# Function to check and print status
check_and_print() {
  if [ $? -eq 0 ]; then
    echo "[OK] $1"
  else
    echo "[FAIL] $1"
    exit 1
  fi
}

# Check config files
[ -f "${CONFIG_DIR}/fit_params_emp_hier_${FIT_CONFIG}.conf" ]
check_and_print "Fit config file: ${CONFIG_DIR}/fit_params_emp_hier_${FIT_CONFIG}.conf"

[ -f "${CONFIG_DIR}/fit_params_emp_indiv_${FIT_CONFIG}.conf" ]
check_and_print "Fit config file: ${CONFIG_DIR}/fit_params_emp_indiv_${FIT_CONFIG}.conf"

[ -f "${CONFIG_DIR}/data_params_${DATA_CONFIG}.conf" ]
check_and_print "Data config file: ${CONFIG_DIR}/data_params_${DATA_CONFIG}.conf"

# Check R script
[ -f "$R_SCRIPT" ]
check_and_print "R script: $R_SCRIPT"

# Check R script 2
[ -f "$R_SCRIPT2" ]
check_and_print "R script 2: $R_SCRIPT2"

# Check R script 3
[ -f "$R_SCRIPT3" ]
check_and_print "R script 3: $R_SCRIPT3"

# Check SBatch script
[ -f "$SBATCH_SCRIPT" ]
check_and_print "SBATCH script: $SBATCH_SCRIPT"

# Check data file
check_and_print "Model type: $MODEL_TYPE"

# Check for required R packages
Rscript -e "if (!all(c('rstan', 'optparse') %in% installed.packages()[,'Package'])) quit(status=1)"
check_and_print "Required R packages installed"

# Convert comma-separated model names to array
IFS=',' read -ra MODEL_ARRAY <<< "$MODEL_NAMES"

# Check model files and submit jobs
for MODEL_NAME in "${MODEL_ARRAY[@]}"; do

  if $DRY_RUN; then
    echo "Dry run for model ${TASK}_${GROUP_TYPE}_${MODEL_NAME}_${MODEL_TYPE}:"
    echo "  SLURM job would be submitted with:"
    echo "    MODEL_NAME=$MODEL_NAME"
    echo "    FIT_CONFIG=$FIT_CONFIG"
    echo "    DATA_CONFIG=$DATA_CONFIG"
    echo "    USER_EMAIL=$USER_EMAIL"
    echo "    MODEL_TYPE=$MODEL_TYPE"
    echo "    TASK=$TASK"
    echo "    GROUP_TYPE=$GROUP_TYPE"
    echo "    CHECK_ITER=$CHECK_ITER"
    echo "    SUBS_FILE=$SUBS_FILE"
    
    # Actually run the R script in dry-run mode
    source "${CONFIG_DIR}/fit_params_emp_hier_${FIT_CONFIG}.conf"
    source "${CONFIG_DIR}/data_params_${DATA_CONFIG}.conf"
    source "${CONFIG_DIR}/fit_params_emp_indiv_${FIT_CONFIG}.conf"
  else
    JOB_NAME="FIT-empbayes_task-${TASK}_cohort-${SOURCE}_model-${MODEL_NAME}"
    job_id=$(sbatch --parsable \
      --job-name=$JOB_NAME \
      --export=JOB_NAME=$JOB_NAME,MODEL_NAME=$MODEL_NAME,FIT_CONFIG=$FIT_CONFIG,DATA_CONFIG=$DATA_CONFIG,USER_EMAIL=$USER_EMAIL,MODEL_TYPE=$MODEL_TYPE,TASK=$TASK,GROUP_TYPE=$GROUP_TYPE,CHECK_ITER=$CHECK_ITER,STEPS=$STEPS,SUBS_FILE=$SUBS_FILE,SOURCE=$SOURCE,SES=$SES\
      $SBATCH_SCRIPT)
    if [ $? -eq 0 ]; then
      echo "Submitted job for model $MODEL_NAME with ID: $job_id"
    else
      echo "Failed to submit job for model $MODEL_NAME"
    fi
  fi
  echo
done

if $DRY_RUN; then
  echo "This was a dry run. No jobs were actually submitted."
else
  echo "Use 'squeue -u $USER' to monitor your jobs."
  echo "Use 'squeue -l --me' to monitor your jobs."
fi