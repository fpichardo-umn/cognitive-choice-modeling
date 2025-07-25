#!/bin/bash

# PR Submission Script
# Submits modeling PR jobs for group/hierarchical models
# Simplified version focused on group/hier types with sensible defaults

# Function to print usage
print_usage() {
  echo "Usage: $0 -k <task> -m <model> -s <source> -e <email> [options]"
  echo "Example: $0 -k igt_mod -m ev -s combined -e your@email.edu"
  echo "Example: $0 -k igt_mod -m ev -s combined -e your@email.edu -g group -c fit,pr_genparams"
  echo "Example: $0 -k igt_mod -m ev -s combined -e your@email.edu --n-subjects 150 --n-iter 3000"
  echo ""
  echo "Required Options:"
  echo "  -k    Task name (e.g., igt_mod)"
  echo "  -m    Model name (e.g., ev, orl)"
  echo "  -s    Data source (e.g., combined, ahrb, adb)"
  echo "  -e    Your email address"
  echo ""
  echo "Optional Options:"
  echo "  -g    Group type: 'group' or 'hier' (default: hier)"
  echo "  -S    Session identifier (default: 00)"
  echo "  -c    Components to run (default: all)"
  echo "        Options: all, fit, pr_genparams, pr_simulate, pr_recovery"
  echo "        Can be comma-separated: fit,pr_genparams,pr_simulate"
  echo "  -n    Dry run (show what would be submitted without submitting)"
  echo ""
  echo "PR Parameter Options:"
  echo "  --n-subjects N    Number of subjects (default: 200)"
  echo "  --n-trials N      Number of trials (default: 100)"
  echo "  --n-iter N        MCMC iterations (default: 5000)"
  echo "  --check-iter N    Checkpoint interval (default: 10000)"
  echo "  --indiv           Enable individual fitting (default: disabled for hier)"
  echo ""
  echo "Defaults:"
  echo "  - Group type: hier"
  echo "  - Components: all"
  echo "  - Session: 00"
  echo "  - Number of subjects: 200"
  echo "  - Number of trials: 100"
  echo "  - MCMC iterations: 5000"
  echo "  - No individual fitting (hierarchical only)"
  echo "  - 72 hour time limit, 64GB memory"
  echo ""
  exit 1
}

# Parse command line arguments with defaults
DRY_RUN=false
GROUP_TYPE="hier"
COMPONENTS="all"
SESSION="00"
N_SUBJECTS=200
N_TRIALS=100
N_ITER=5000
CHECK_ITER=10000
NO_INDIV=true

while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
      print_usage
      exit 0
      ;;
    -k) TASK=$2; shift 2 ;;
    -m) MODEL=$2; shift 2 ;;
    -s) SOURCE=$2; shift 2 ;;
    -e) USER_EMAIL=$2; shift 2 ;;
    -g) GROUP_TYPE=$2; shift 2 ;;
    -S) SESSION=$2; shift 2 ;;
    -c) COMPONENTS=$2; shift 2 ;;
    -n) DRY_RUN=true; shift ;;
    --n-subjects) N_SUBJECTS=$2; shift 2 ;;
    --n-trials) N_TRIALS=$2; shift 2 ;;
    --n-iter) N_ITER=$2; shift 2 ;;
    --check-iter) CHECK_ITER=$2; shift 2 ;;
    --indiv) NO_INDIV=false; shift ;;
    -*) echo "Invalid option: $1" >&2; print_usage ;;
    *) echo "Unexpected argument: $1" >&2; print_usage ;;
  esac
done

# Check required arguments
if [ -z "$TASK" ] || [ -z "$MODEL" ] || [ -z "$SOURCE" ] || [ -z "$USER_EMAIL" ]; then
  echo "Error: Task (-k), model (-m), source (-s), and email (-e) are required."
  print_usage
fi

# Validate group type
if [ "$GROUP_TYPE" != "group" ] && [ "$GROUP_TYPE" != "hier" ]; then
  echo "Error: Group type must be 'group' or 'hier'"
  print_usage
fi

# Source helper bash script for directory functions
SUBMIT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
SCRIPT_DIR="${SUBMIT_DIR}/.."
source "${SCRIPT_DIR}/helpers/dir_helpers.sh"

# Get task-specific directories
PROJ_DIR="${SCRIPT_DIR}/.."
LOG_DIR="$(get_log_dir "parameter_recovery")"
BATCH_SCRIPT="${SCRIPT_DIR}/sbatch_resc_scripts/resources_pr.sbatch"

# Generate job name with BIDS-style format
JOB_NAME="PR_task-${TASK}_model-${MODEL}_group-${GROUP_TYPE}_comp-${COMPONENTS//,/-}"
if [ ! -z "$SESSION" ]; then
  JOB_NAME="${JOB_NAME}_ses-${SESSION}"
fi

# Build environment variables for the batch script (use semicolon to avoid comma conflicts)
COMPONENTS_ENCODED=$(echo "$COMPONENTS" | sed 's/,/;/g')
ENV_VARS="ALL,TASK=$TASK,MODEL=$MODEL,SOURCE=$SOURCE,GROUP_TYPE=$GROUP_TYPE"
ENV_VARS="${ENV_VARS},COMPONENTS=$COMPONENTS_ENCODED"
ENV_VARS="${ENV_VARS},SESSION=$SESSION,N_SUBJECTS=$N_SUBJECTS,N_TRIALS=$N_TRIALS,N_ITER=$N_ITER,CHECK_ITER=$CHECK_ITER,NO_INDIV=$NO_INDIV"

if $DRY_RUN; then
  echo "====== DRY RUN: WOULD SUBMIT WITH THESE PARAMETERS ======"
  echo "Task: $TASK"
  echo "Model: $MODEL"
  echo "Source: $SOURCE"
  echo "Group type: $GROUP_TYPE"
  echo "Components: $COMPONENTS"
  echo "Session: $SESSION"
  echo "Number of subjects: $N_SUBJECTS"
  echo "Number of trials: $N_TRIALS"
  echo "MCMC iterations: $N_ITER"
  echo "Checkpoint interval: $CHECK_ITER"
  echo "No individual fitting: $NO_INDIV"
  echo "Email: $USER_EMAIL"
  echo "Job name: $JOB_NAME"
  echo "Log directory: $LOG_DIR"
  echo ""
  echo "Command that would be executed:"
  echo "sbatch --parsable \\"
  echo "  --job-name=$JOB_NAME \\"
  echo "  --mail-user=$USER_EMAIL \\"
  echo "  --export=$ENV_VARS \\"
  echo "  --output=\"${LOG_DIR}/${JOB_NAME}_%j.out\" \\"
  echo "  --error=\"${LOG_DIR}/${JOB_NAME}_%j.err\" \\"
  echo "  $BATCH_SCRIPT"
  echo ""
  echo "PR command that would be run:"
  echo "run_full_PR.sh -k $TASK -m $MODEL -s $SOURCE -g $GROUP_TYPE --ses $SESSION \\"
  echo "  --n-subjects $N_SUBJECTS --n-trials $N_TRIALS --n-iter $N_ITER --check-iter $CHECK_ITER \\"
  if [ "$NO_INDIV" = true ]; then
    echo "  --no-indiv $COMPONENTS"
  else
    echo "  $COMPONENTS"
  fi
  echo "====== DRY RUN COMPLETED ======"
else
  # Ensure the necessary directories exist
  ensure_dir_exists "$LOG_DIR"
  
  # Submit the job
  job_id=$(sbatch --parsable \
    --job-name=$JOB_NAME \
    --mail-user=$USER_EMAIL \
    --export=$ENV_VARS \
    --output="${LOG_DIR}/${JOB_NAME}_%j.out" \
    --error="${LOG_DIR}/${JOB_NAME}_%j.err" \
    $BATCH_SCRIPT)
  
  if [ $? -eq 0 ]; then
    echo "Submitted PR job with ID: $job_id"
    echo "Job name: $JOB_NAME"
    echo "Task: $TASK, Model: $MODEL, Source: $SOURCE"
    echo "Group type: $GROUP_TYPE, Components: $COMPONENTS"
    echo "Session: $SESSION, Subjects: $N_SUBJECTS, Trials: $N_TRIALS"
    echo "MCMC iterations: $N_ITER, Checkpoint: $CHECK_ITER"
    echo "Check logs in: ${LOG_DIR}/"
  else
    echo "Failed to submit PR job"
    exit 1
  fi
fi