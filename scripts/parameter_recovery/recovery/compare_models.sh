#!/bin/bash

# compare_models.sh - Script to run model comparison for mIGT models
# Usage: ./compare_models.sh [task] [models] [group] [cohort] [session]

# Default values
TASK=${1:-"igt"}
MODELS=${2:-"ev,pvl,baseline"}  # Comma-separated list of models
GROUP=${3:-"sing"}  # sing, hier, batch
COHORT=${4:-"sim"}
SESSION=${5:-""}

# Create directory for the script if it doesn't exist
SCRIPT_DIR="$(dirname "$0")"
mkdir -p "$SCRIPT_DIR"

# Print parameters
echo "Running model comparison with parameters:"
echo "  Task: $TASK"
echo "  Models: $MODELS"
echo "  Group: $GROUP"
echo "  Cohort: $COHORT"
echo "  Session: $SESSION"

# Run model comparison
Rscript "$SCRIPT_DIR/compare_models.R" \
  --task "$TASK" \
  --models "$MODELS" \
  --group "$GROUP" \
  --cohort "$COHORT" \
  --session "$SESSION" \
  --render

echo "Model comparison complete."
