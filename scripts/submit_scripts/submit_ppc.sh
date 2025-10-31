#!/bin/bash

# PPC Submission Script
# Submits PPC jobs with config file support

# Function to print usage
print_usage() {
  echo "Usage: $0 -m <model_name> -k <task> -s <cohort> -e <email> [options]"
  echo ""
  echo "Example: $0 -m ev -k igt_mod -s combined -g batch_001 -e your@email.edu"
  echo "Example: $0 -m ev -k igt_mod -s combined -g hier -e your@email.edu -F complex -D igt_mod"
  echo ""
  echo "Required Options:"
  echo "  -m    Model name (e.g., ev, orl, ev_ddm)"
  echo "  -k    Task name (e.g., igt_mod)"
  echo "  -s    Source/cohort identifier (e.g., combined, ahrb)"
  echo "  -e    Your email address"
  echo ""
  echo "Optional Options:"
  echo "  -g    Group name/ID (e.g., batch_001, hier) (default: batch_001)"
  echo "  -S    Session identifier (default: 00)"
  echo "  -F    Fit config name (default: default)"
  echo "  -D    Data config name (default: default)"
  echo "  -n    Number of simulations per subject (default: 100)"
  echo "  -b    Block size (default: 20)"
  echo "  -p    Sampling method (default: weighted, options: random, width, weighted)"
  echo "  -w    Width control for width sampling (default: 0.95)"
  echo "  -x    Excludes file path (default: Data/txt/subs/subject_ids_excludes.txt)"
  echo "  -i    Information criteria method (default: loo, options: loo, waic)"
  echo "  -t    Steps to run (default: all, options: comma-separated list of simulate,stats,loglik,report,all)"
  echo "  -d    Dry run (optional)"
  echo ""
  echo "Note: RT parameters (rt_method, RTbound_min_ms, RTbound_max_ms) are sourced from data config file"
  echo ""
  exit 1
}

# --- Set Defaults ---
DRY_RUN=false
GROUP_NAME="batch_001"
SESSION="00"
FIT_CONFIG="default"
DATA_CONFIG="default"
N_SIMS=100
BLOCK_SIZE=20
SAMPLING="weighted"
WIDTH_CONTROL=0.95
EXCLUDES_FILE=""
IC_METHOD="loo"
STEPS="all"

# --- Parse Command Line Arguments ---
while getopts ":m:k:s:e:g:S:F:D:n:b:p:w:x:i:t:dh" opt; do
  case $opt in
    m) MODEL_NAME=$OPTARG ;;
    k) TASK=$OPTARG ;;
    s) COHORT=$OPTARG ;;
    e) USER_EMAIL=$OPTARG ;;
    g) GROUP_NAME=$OPTARG ;;
    S) SESSION=$OPTARG ;;
    F) FIT_CONFIG=$OPTARG ;;
    D) DATA_CONFIG=$OPTARG ;;
    n) N_SIMS=$OPTARG ;;
    b) BLOCK_SIZE=$OPTARG ;;
    p) SAMPLING=$OPTARG ;;
    w) WIDTH_CONTROL=$OPTARG ;;
    x) EXCLUDES_FILE=$OPTARG ;;
    i) IC_METHOD=$OPTARG ;;
    t) STEPS=$OPTARG ;;
    d) DRY_RUN=true ;;
    h) print_usage ;;
    \?) echo "Invalid option -$OPTARG" >&2; print_usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; print_usage ;;
  esac
done

# --- Validate Required Arguments ---
if [ -z "$MODEL_NAME" ] || [ -z "$TASK" ] || [ -z "$COHORT" ] || [ -z "$USER_EMAIL" ]; then
  echo "Error: Model name (-m), task (-k), cohort (-s), and email (-e) are required."
  print_usage
fi

# --- Set up Directories and Paths ---
SUBMIT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
SCRIPT_DIR="${SUBMIT_DIR}/.."
PROJ_DIR="${SCRIPT_DIR}/.."

# Source directory helper if available
if [ -f "${SCRIPT_DIR}/helpers/dir_helpers.sh" ]; then
  source "${SCRIPT_DIR}/helpers/dir_helpers.sh"
  LOG_DIR="$(get_log_dir "ppc")"
else
  LOG_DIR="${PROJ_DIR}/log_files/ppc"
fi

BATCH_SCRIPT="${SCRIPT_DIR}/sbatch_resc_scripts/resources_ppc.sbatch"

# --- Generate Job Name ---
JOB_NAME="PPC_task-${TASK}_cohort-${COHORT}_model-${MODEL_NAME}_group-${GROUP_NAME}"
if [ ! -z "$SESSION" ] && [ "$SESSION" != "00" ]; then
  JOB_NAME="${JOB_NAME}_ses-${SESSION}"
fi

# --- Dry Run or Submit ---
if $DRY_RUN; then
  echo "====== DRY RUN: WOULD SUBMIT WITH THESE PARAMETERS ======"
  echo "Job name: $JOB_NAME"
  echo "Email: $USER_EMAIL"
  echo "Log directory: $LOG_DIR"
  echo "---"
  echo "Task: $TASK"
  echo "Model: $MODEL_NAME"
  echo "Cohort: $COHORT"
  echo "Group name: $GROUP_NAME"
  echo "Session: $SESSION"
  echo "---"
  echo "Config files:"
  echo "  Fit config: $FIT_CONFIG"
  echo "  Data config: $DATA_CONFIG"
  echo "---"
  echo "PPC parameters:"
  echo "  Number of simulations: $N_SIMS"
  echo "  Block size: $BLOCK_SIZE"
  echo "  Sampling method: $SAMPLING"
  echo "  Width control: $WIDTH_CONTROL"
  echo "  Excludes file: $EXCLUDES_FILE"
  echo "  IC method: $IC_METHOD"
  echo "  Steps to run: $STEPS"
  echo ""
  echo "Note: RT parameters will be sourced from: scripts/configs/data_params_${DATA_CONFIG}.conf"
  echo ""
  echo "Command that would be executed:"
  echo "sbatch --parsable \\"
  echo "  --job-name=$JOB_NAME \\"
  echo "  --mail-user=$USER_EMAIL \\"
  echo "  --output=\"${LOG_DIR}/${JOB_NAME}_%j.out\" \\"
  echo "  --error=\"${LOG_DIR}/${JOB_NAME}_%j.err\" \\"
  echo "  --export=ALL,MODEL_NAME=$MODEL_NAME,TASK=$TASK,COHORT=$COHORT,GROUP_NAME=$GROUP_NAME,SESSION=$SESSION,FIT_CONFIG=$FIT_CONFIG,DATA_CONFIG=$DATA_CONFIG,N_SIMS=$N_SIMS,BLOCK_SIZE=$BLOCK_SIZE,SAMPLING=$SAMPLING,WIDTH_CONTROL=$WIDTH_CONTROL,EXCLUDES_FILE=$EXCLUDES_FILE,IC_METHOD=$IC_METHOD,STEPS=$STEPS \\"
  echo "  $BATCH_SCRIPT"
  echo "====== DRY RUN COMPLETED ======"
else
  # Ensure the necessary directories exist
  mkdir -p "$LOG_DIR"

  # Submit the job
  job_id=$(sbatch --parsable \
    --job-name=$JOB_NAME \
    --mail-user=$USER_EMAIL \
    --output="${LOG_DIR}/${JOB_NAME}_%j.out" \
    --error="${LOG_DIR}/${JOB_NAME}_%j.err" \
    --export=ALL,MODEL_NAME=$MODEL_NAME,TASK=$TASK,COHORT=$COHORT,GROUP_NAME=$GROUP_NAME,SESSION=$SESSION,FIT_CONFIG=$FIT_CONFIG,DATA_CONFIG=$DATA_CONFIG,N_SIMS=$N_SIMS,BLOCK_SIZE=$BLOCK_SIZE,SAMPLING=$SAMPLING,WIDTH_CONTROL=$WIDTH_CONTROL,EXCLUDES_FILE=$EXCLUDES_FILE,IC_METHOD=$IC_METHOD,STEPS=$STEPS \
    $BATCH_SCRIPT)

  if [ $? -eq 0 ]; then
    echo "Submitted PPC job with ID: $job_id"
    echo "   Job name: $JOB_NAME"
    echo "   Task: $TASK, Model: $MODEL_NAME, Cohort: $COHORT"
    echo "   Group: $GROUP_NAME, Session: $SESSION"
    echo "   Fit config: $FIT_CONFIG, Data config: $DATA_CONFIG"
    echo "   Check logs in: ${LOG_DIR}/"
  else
    echo "Failed to submit PPC job"
    exit 1
  fi
fi
