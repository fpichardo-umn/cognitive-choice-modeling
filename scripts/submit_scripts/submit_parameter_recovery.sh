#!/bin/bash

# Function to print usage
print_usage() {
  echo "Usage: $0 -m <model_name> -k <task> -g <batch_name> -e <email> [-n <n_subjects>] [-d <data_config>] [-s <seed>] [-x <excludes_file>] [-i <n_iterations>] [-p <n_ppc>] [-D] [-r]"
  echo "Example: $0 -m pvlbothsd -k igt_mod -g batch_001 -e your@email.edu"
  echo "Options:"
  echo "  -m    Model name (e.g., pvlbothsd)"
  echo "  -k    Task name (e.g., igt_mod)"
  echo "  -g    Batch name/ID (e.g., batch_001)"
  echo "  -e    Your email address"
  echo "  -n    Number of subjects to simulate (default: 200)"
  echo "  -d    Data parameters config name (default: mbSPSepse)"
  echo "  -s    Random seed (default: 12345)"
  echo "  -x    Excludes file path (default: Data/txt/subs/subject_ids_excludes.txt)"
  echo "  -i    Number of iterations for PPC (default: 1000)"
  echo "  -p    Number of samples for PPC (default: 100)"
  echo "  -D    Dry run (optional)"
  echo "  -r    Run PPC only (skips parameter generation, simulation, and recovery)"
  exit 1
}

# Parse command line arguments
DRY_RUN=false
PPC_ONLY=false
N_SUBJECTS=200
DATA_CONFIG="mbSPSepse"
SEED=12345
EXCLUDES_FILE="Data/txt/subs/subject_ids_excludes.txt"
N_ITERATIONS=1000
N_PPC=100

while getopts ":m:k:g:e:n:d:s:x:i:p:Dr" opt; do
  case $opt in
    m) MODEL_NAME=$OPTARG ;;
    k) TASK=$OPTARG ;;
    g) BATCH_NAME=$OPTARG ;;
    e) USER_EMAIL=$OPTARG ;;
    n) N_SUBJECTS=$OPTARG ;;
    d) DATA_CONFIG=$OPTARG ;;
    s) SEED=$OPTARG ;;
    x) EXCLUDES_FILE=$OPTARG ;;
    i) N_ITERATIONS=$OPTARG ;;
    p) N_PPC=$OPTARG ;;
    D) DRY_RUN=true ;;
    r) PPC_ONLY=true ;;
    \?) echo "Invalid option -$OPTARG" >&2; print_usage ;;
  esac
done

# Check required arguments
if [ -z "$MODEL_NAME" ] || [ -z "$TASK" ] || [ -z "$BATCH_NAME" ] || [ -z "$USER_EMAIL" ]; then
  echo "Error: Model name, task, batch name, and email are required."
  print_usage
fi

# Set up directories
SUBMIT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
SCRIPT_DIR="${SUBMIT_DIR}/.."
PROJ_DIR="${SCRIPT_DIR}/.."
BATCH_SCRIPT="${SCRIPT_DIR}/sbatch_resc_scripts/resources_parameter_recovery.sbatch"

# Generate parameter file path
PARAMS_FILE="Data/sim/params/${TASK}_${BATCH_NAME}_${MODEL_NAME}_subj-${N_SUBJECTS}_desc-sim_params_${DATA_CONFIG}.rds"
SIM_DATA_FILE="Data/sim/rds/${TASK}_${BATCH_NAME}_${MODEL_NAME}_desc-sim_params.rds"
FIT_OUTPUT_FILE="Data/sim/fit/fits_${TASK}_${MODEL_NAME}_${BATCH_NAME}.rds"

if $DRY_RUN; then
  echo "Dry run - would submit with these parameters:"
  echo "Model name: $MODEL_NAME"
  echo "Task: $TASK"
  echo "Batch name: $BATCH_NAME"
  echo "Number of subjects: $N_SUBJECTS"
  echo "Data config: $DATA_CONFIG"
  echo "Seed: $SEED"
  echo "Excludes file: $EXCLUDES_FILE"
  echo "PPC iterations: $N_ITERATIONS"
  echo "PPC samples: $N_PPC"
  echo "Email: $USER_EMAIL"
  echo "PPC only mode: $PPC_ONLY"
else
  JOB_NAME="pr_${TASK}_${MODEL_NAME}_${BATCH_NAME}"
  job_id=$(sbatch --parsable \
    --job-name=$JOB_NAME \
    --mail-user=$USER_EMAIL \
    --export=ALL,MODEL_NAME=$MODEL_NAME,TASK=$TASK,BATCH_NAME=$BATCH_NAME,N_SUBJECTS=$N_SUBJECTS,DATA_CONFIG=$DATA_CONFIG,SEED=$SEED,EXCLUDES_FILE=$EXCLUDES_FILE,N_ITERATIONS=$N_ITERATIONS,N_PPC=$N_PPC,PPC_ONLY=$PPC_ONLY,PARAMS_FILE=$PARAMS_FILE,SIM_DATA_FILE=$SIM_DATA_FILE,FIT_OUTPUT_FILE=$FIT_OUTPUT_FILE \
    $BATCH_SCRIPT)
  
  if [ $? -eq 0 ]; then
    echo "Submitted parameter recovery job with ID: $job_id"
    echo "Model: $MODEL_NAME, Task: $TASK, Batch: $BATCH_NAME"
  else
    echo "Failed to submit parameter recovery job"
  fi
fi