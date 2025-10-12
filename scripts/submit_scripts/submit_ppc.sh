#!/bin/bash

# Function to print usage
print_usage() {
  echo "Usage: $0 -m <model_name> -k <task> -g <group_name> -e <email> [options]"
  echo "Example: $0 -m ev_ddm -k igt_mod -g batch_001 -e your@email.edu"
  echo "Options:"
  echo "  -m    Model name (e.g., ev_ddm, pvlbothsd)"
  echo "  -k    Task name (e.g., igt_mod)"
  echo "  -g    Group name/ID (e.g., batch_001)"
  echo "  -e    Your email address"
  echo "  -n    Number of simulations per subject (default: 100)"
  echo "  -b    Block size (default: 20)"
  echo "  -s    Sampling method (default: weighted, options: random, width, weighted)"
  echo "  -w    Width control for width sampling (default: 0.95)"
  echo "  -r    RT method (default: remove, options: all, remove, force, adaptive)"
  echo "  -R    RT lower bound in milliseconds (default: 100)"
  echo "  -U    RT upper bound in milliseconds (default: 2500)"
  echo "  -x    Excludes file path (default: Data/txt/subs/subject_ids_excludes.txt)"
  echo "  -i    Information criteria method (default: loo, options: loo, waic)"
  echo "  -t    Steps to run (default: all, options: comma-separated list of simulate,stats,loglik,report,all)"
  echo "  -D    Dry run (optional)"
  exit 1
}

# Parse command line arguments
DRY_RUN=false
N_SIMS=100
BLOCK_SIZE=20
SAMPLING="weighted"
WIDTH_CONTROL=0.95
RT_METHOD="remove"
RT_BOUND_MIN_MS=100
RT_BOUND_MAX_MS=2500
EXCLUDES_FILE="Data/txt/subs/subject_ids_excludes.txt"
IC_METHOD="loo"
STEPS="all"

while getopts ":m:k:g:e:n:b:s:w:r:R:U:x:i:t:D" opt; do
  case $opt in
    m) MODEL_NAME=$OPTARG ;;
    k) TASK=$OPTARG ;;
    g) GROUP_NAME=$OPTARG ;;
    e) USER_EMAIL=$OPTARG ;;
    n) N_SIMS=$OPTARG ;;
    b) BLOCK_SIZE=$OPTARG ;;
    s) SAMPLING=$OPTARG ;;
    w) WIDTH_CONTROL=$OPTARG ;;
    r) RT_METHOD=$OPTARG ;;
    R) RT_BOUND_MIN_MS=$OPTARG ;;
    U) RT_BOUND_MAX_MS=$OPTARG ;;
    x) EXCLUDES_FILE=$OPTARG ;;
    i) IC_METHOD=$OPTARG ;;
    t) STEPS=$OPTARG ;;
    D) DRY_RUN=true ;;
    \?) echo "Invalid option -$OPTARG" >&2; print_usage ;;
  esac
done

# Check required arguments
if [ -z "$MODEL_NAME" ] || [ -z "$TASK" ] || [ -z "$GROUP_NAME" ] || [ -z "$USER_EMAIL" ]; then
  echo "Error: Model name, task, group name, and email are required."
  print_usage
fi

# Set up directories
SUBMIT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
SCRIPT_DIR="${SUBMIT_DIR}/.."
PROJ_DIR="${SCRIPT_DIR}/.."
BATCH_SCRIPT="${SCRIPT_DIR}/sbatch_resc_scripts/resources_ppc.sbatch"

# Define output file paths
SIM_FILE="Data/ppc/simulations/ppc_simulations_task-${TASK}_model-${MODEL_NAME}_group-${GROUP_NAME}.rds"
STATS_FILE="Data/ppc/stats/ppc_summary_task-${TASK}_model-${MODEL_NAME}_group-${GROUP_NAME}.csv"
LOGLIK_FILE="Data/ppc/loglik/ppc_loglik_task-${TASK}_model-${MODEL_NAME}_group-${GROUP_NAME}.rds"
REPORT_FILE="Data/ppc/reports/ppc_report_task-${TASK}_model-${MODEL_NAME}_group-${GROUP_NAME}.html"

if $DRY_RUN; then
  echo "Dry run - would submit with these parameters:"
  echo "Model name: $MODEL_NAME"
  echo "Task: $TASK"
  echo "Group name: $GROUP_NAME"
  echo "Number of simulations: $N_SIMS"
  echo "Block size: $BLOCK_SIZE"
  echo "Sampling method: $SAMPLING"
  echo "Width control: $WIDTH_CONTROL"
  echo "RT method: $RT_METHOD"
  echo "RT lower bound (ms): $RT_BOUND_MIN_MS"
  echo "RT upper bound (ms): $RT_BOUND_MAX_MS"
  echo "Excludes file: $EXCLUDES_FILE"
  echo "IC method: $IC_METHOD"
  echo "Steps to run: $STEPS"
  echo "Email: $USER_EMAIL"
else
  JOB_NAME="ppc_${TASK}_${MODEL_NAME}_${GROUP_NAME}"
  job_id=$(sbatch --parsable \
    --job-name=$JOB_NAME \
    --mail-user=$USER_EMAIL \
    --export=ALL,MODEL_NAME=$MODEL_NAME,TASK=$TASK,GROUP_NAME=$GROUP_NAME,N_SIMS=$N_SIMS,BLOCK_SIZE=$BLOCK_SIZE,SAMPLING=$SAMPLING,WIDTH_CONTROL=$WIDTH_CONTROL,RT_METHOD=$RT_METHOD,RT_BOUND_MIN_MS=$RT_BOUND_MIN_MS,RT_BOUND_MAX_MS=$RT_BOUND_MAX_MS,EXCLUDES_FILE=$EXCLUDES_FILE,IC_METHOD=$IC_METHOD,STEPS=$STEPS \
    $BATCH_SCRIPT)
  
  if [ $? -eq 0 ]; then
    echo "Submitted PPC job with ID: $job_id"
    echo "Model: $MODEL_NAME, Task: $TASK, Group: $GROUP_NAME"
    echo "Expected output files:"
    echo "  Simulations: $SIM_FILE"
    echo "  Statistics: $STATS_FILE"
    echo "  Log-likelihood: $LOGLIK_FILE"
    echo "  Report: $REPORT_FILE"
  else
    echo "Failed to submit PPC job"
  fi
fi