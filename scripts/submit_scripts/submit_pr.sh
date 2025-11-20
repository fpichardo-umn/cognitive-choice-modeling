#!/bin/bash

# PR Submission Script
# Submits modeling PR jobs for group/hierarchical models
# Refactored to align with the structure of submit_ppc.sh

# Function to print usage
print_usage() {
  echo "Usage: $0 -k <task> -m <model> -s <source> -e <email> [options]"
  echo ""
  echo "Example: $0 -k igt_mod -m ev -s combined -e your@email.edu -g group -c fit,pr_genparams"
  echo "Example: $0 -k igt_mod -m ev -s combined -e your@email.edu -N 150 -F emp_hier_complex"
  echo ""
  echo "Required Options:"
  echo "  -k    Task name (e.g., igt_mod)"
  echo "  -m    Model name (e.g., ev, orl)"
  echo "  -s    Data source (e.g., combined, ahrb, adb)"
  echo "  -e    Your email address"
  echo ""
  echo "Optional Options:"
  echo "  -g    Group type: 'sing' (for batch), 'group' or 'hier' (default: hier)"
  echo "  -S    Session identifier (default: 00)"
  echo "  -c    Components to run, comma-separated (default: all)"
  echo "        (options: fit, pr_genparams, pr_simulate, pr_recovery, all)"
  echo "  -i    Enable individual fitting (default: disabled for hier)"
  echo "  -d    Dry run (show what would be submitted without submitting)"
  echo "  -n    Number of subjects to fit (default: 200)"
  echo ""
  echo "Configuration & PR Options:"
  echo "  -F    Fit config name (default: default)"
  echo "  -I    Sim fit config name for recovery (default: sim)"
  echo "  -D    Data config name (default: default)"
  echo "  -M    Parameter generation method (default: ibSPSepse)"
  echo "  -f    Subject file in Data/raw/{source}/ses-{ses}/ (default: subject_ids_all.txt)"
  echo "  -N    Number of subjects for PR (default: 200)"
  echo "  -T    Number of trials for PR (default: 100)"
  echo ""
  echo "Defaults:"
  echo "  - Group type: hier, Session: 00, Components: all"
  echo "  - Subjects: 200, Trials: 100"
  echo "  - Fit config: default, Data config: default, Method: ibSPSepse"
  echo "  - 96 hour time limit, 64GB memory"
  echo ""
  exit 1
}

# --- Set Defaults ---
DRY_RUN=false
GROUP_TYPE="hier"
COMPONENTS="all"
SESSION="00"
N_SUBJECTS_FIT=200
N_SUBJECTS_PR=200
N_TRIALS=100
NO_INDIV=true
FIT_CONFIG="default"
SIM_CONFIG="sim"
DATA_CONFIG="default"
METHOD="ibSPSepse"
SUBS_FILE="subject_ids_all.txt"

# --- Parse Command Line Arguments ---
while getopts ":k:m:s:e:g:S:c:F:I:D:M:f:N:n:T:idh" opt; do
  case $opt in
    k) TASK=$OPTARG ;;
    m) MODEL=$OPTARG ;;
    s) SOURCE=$OPTARG ;;
    e) USER_EMAIL=$OPTARG ;;
    g) GROUP_TYPE=$OPTARG ;;
    S) SESSION=$OPTARG ;;
    c) COMPONENTS=$OPTARG ;;
    F) FIT_CONFIG=$OPTARG ;;
    I) SIM_CONFIG=$OPTARG ;;
    D) DATA_CONFIG=$OPTARG ;;
    M) METHOD=$OPTARG ;;
    f) SUBS_FILE=$OPTARG ;;
    N) N_SUBJECTS_PR=$OPTARG ;;
    n) N_SUBJECTS_FIT=$OPTARG ;;
    T) N_TRIALS=$OPTARG ;;
    i) NO_INDIV=false ;;
    d) DRY_RUN=true ;;
    h) print_usage ;;
    \?) echo "Invalid option -$OPTARG" >&2; print_usage ;;
    :) echo "Option -$OPTARG requires an argument." >&2; print_usage ;;
  esac
done

# --- Validate Required Arguments ---
if [ -z "$TASK" ] || [ -z "$MODEL" ] || [ -z "$SOURCE" ] || [ -z "$USER_EMAIL" ]; then
  echo "Error: Task (-k), model (-m), source (-s), and email (-e) are required."
  print_usage
fi

if [ "$GROUP_TYPE" != "sing" ] && [ "$GROUP_TYPE" != "group" ] && [ "$GROUP_TYPE" != "hier" ]; then
  echo "Error: Group type (-g) must be 'sing' (for batch), 'group' or 'hier'"
  print_usage
fi

# For batch mode (sing), individual fitting is required
if [ "$GROUP_TYPE" == "sing" ]; then
   NO_INDIV=false
fi

# --- Set up Directories and Paths ---
SUBMIT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
SCRIPT_DIR="${SUBMIT_DIR}/.."
source "${SCRIPT_DIR}/helpers/dir_helpers.sh" # Assumes this helper exists and works

LOG_DIR="$(get_log_dir "parameter_recovery")"
BATCH_SCRIPT="${SCRIPT_DIR}/sbatch_resc_scripts/resources_pr.sbatch"

# --- Generate Job Name ---
JOB_NAME="PR_task-${TASK}_cohort-${SOURCE}_model-${MODEL}_group-${GROUP_TYPE}_comp-${COMPONENTS//,/-}"
if [ ! -z "$SESSION" ]; then
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
  echo "Model: $MODEL"
  echo "Source: $SOURCE"
  echo "Group type: $GROUP_TYPE"
  echo "Components: $COMPONENTS"
  echo "Session: $SESSION"
  echo "---"
  echo "Fit config: $FIT_CONFIG"
  echo "Sim fit config: $SIM_CONFIG"
  echo "Data config: $DATA_CONFIG"
  echo "Fit Subjects: $N_SUBJECTS_FIT"
  echo "PR Method: $METHOD"
  echo "PR Subjects: $N_SUBJECTS_PR"
  echo "PR Trials: $N_TRIALS"
  echo "Individual fitting: $(if [ "$NO_INDIV" = false ]; then echo "enabled"; else echo "disabled"; fi)"
  echo "Subject file: ${SUBS_FILE:-'(not specified - will use default: subject_ids_all.txt)'}"
  echo ""
  echo "Command that would be executed:"
  echo "sbatch --parsable \\"
  echo "  --job-name=$JOB_NAME \\"
  echo "  --mail-user=$USER_EMAIL \\"
  echo "  --output=\"${LOG_DIR}/${JOB_NAME}_%j.out\" \\"
  echo "  --error=\"${LOG_DIR}/${JOB_NAME}_%j.err\" \\"
  echo "  --export=ALL,TASK=$TASK,MODEL=$MODEL,SOURCE=$SOURCE,USER_EMAIL=$USER_EMAIL,GROUP_TYPE=$GROUP_TYPE,COMPONENTS=${COMPONENTS//,/;},SESSION=$SESSION,N_SUBJECTS_FIT=$N_SUBJECTS_FIT,N_SUBJECTS_PR=$N_SUBJECTS_PR,N_TRIALS=$N_TRIALS,NO_INDIV=$NO_INDIV,FIT_CONFIG=$FIT_CONFIG,SIM_CONFIG=$SIM_CONFIG,DATA_CONFIG=$DATA_CONFIG,METHOD=$METHOD,SUBS_FILE=$SUBS_FILE \\"
  echo "  $BATCH_SCRIPT"
  echo "====== DRY RUN COMPLETED ======"
else
  # Ensure the necessary directories exist
  ensure_dir_exists "$LOG_DIR"

  # Submit the job
  job_id=$(sbatch --parsable \
    --job-name=$JOB_NAME \
    --mail-user=$USER_EMAIL \
    --output="${LOG_DIR}/${JOB_NAME}_%j.out" \
    --error="${LOG_DIR}/${JOB_NAME}_%j.err" \
    --export=ALL,TASK=$TASK,MODEL=$MODEL,SOURCE=$SOURCE,USER_EMAIL=$USER_EMAIL,GROUP_TYPE=$GROUP_TYPE,COMPONENTS=${COMPONENTS//,/;},SESSION=$SESSION,N_SUBJECTS_FIT=$N_SUBJECTS_FIT,N_SUBJECTS_PR=$N_SUBJECTS_PR,N_TRIALS=$N_TRIALS,NO_INDIV=$NO_INDIV,FIT_CONFIG=$FIT_CONFIG,SIM_CONFIG=$SIM_CONFIG,DATA_CONFIG=$DATA_CONFIG,METHOD=$METHOD,SUBS_FILE=$SUBS_FILE \
    $BATCH_SCRIPT)

  if [ $? -eq 0 ]; then
    echo "Submitted PR job with ID: $job_id"
    echo "   Job name: $JOB_NAME"
    echo "   Task: $TASK, Model: $MODEL, Source: $SOURCE"
    echo "   Group type: $GROUP_TYPE, Components: $COMPONENTS"
    echo "   Fit Subjects: $N_SUBJECTS_FIT, PR Subjects: $N_SUBJECTS_PR, Trials: $N_TRIALS"
    echo "   Check logs in: ${LOG_DIR}/"
  else
    echo "Failed to submit PR job"
    exit 1
  fi
fi