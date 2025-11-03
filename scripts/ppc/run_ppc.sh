#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 --model MODEL --cohort COHORT [options]

Required arguments:
  --model MODEL          Model name (e.g., ev, pvldelta)
  --cohort COHORT        Cohort identifier

Optional arguments:
  --task TASK            Task name (default: igt_mod)
  --ses SES              Session identifier
  --group GROUP          Fit type: sing (individual) or hier (hierarchical) (default: sing)
  --group_name NAME      Batch identifier (default: batch_001)
  --fit_file FILE        Path to fit file
  --n_sims N             Number of simulations (default: 100)
  --block_size N         Number of trials per block (default: 20)
  --output_dir DIR       Output directory
  --parallel             Use parallel processing
  --n_cores N            Number of cores (default: 2)
  --steps STEPS          Steps to run (default: all)
  --force                Force re-running steps
  --render               Render the Rmd file to HTML
  --sampling METHOD      Sampling method (default: random)
  --width_control NUM    Width control parameter (default: 0.95)
  --rt_method METHOD     RT handling method (default: remove)
  --RTbound_min_ms N     RT lower bound ms (default: 100)
  --RTbound_max_ms N     RT upper bound ms (default: 2500)
  --exclude_file FILE    Path to exclude file
  --ic_method METHOD     Information criterion method (default: loo)
  --dry-run              Show the command but do not execute
  -h, --help             Show this help message
EOF
}

# Defaults
task="igt_mod"
group="sing"
group_name="batch_001"
n_sims=100
block_size=20
parallel=false
n_cores=2
steps="all"
force=false
render=false
sampling="random"
width_control=0.95
rt_method="remove"
RTbound_min_ms=100
RTbound_max_ms=2500
fit_file=""
output_dir=""
ses=""
exclude_file=""
ic_method="loo"
dry_run=false

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --model) model="$2"; shift 2 ;;
    --cohort) cohort="$2"; shift 2 ;;
    --task) task="$2"; shift 2 ;;
    --ses) ses="$2"; shift 2 ;;
    --group) group="$2"; shift 2 ;;
    --group_name) group_name="$2"; shift 2 ;;
    --fit_file) fit_file="$2"; shift 2 ;;
    --n_sims) n_sims="$2"; shift 2 ;;
    --block_size) block_size="$2"; shift 2 ;;
    --output_dir) output_dir="$2"; shift 2 ;;
    --parallel) parallel=true; shift ;;
    --n_cores) n_cores="$2"; shift 2 ;;
    --steps) steps="$2"; shift 2 ;;
    --force) force=true; shift ;;
    --render) render=true; shift ;;
    --sampling) sampling="$2"; shift 2 ;;
    --width_control) width_control="$2"; shift 2 ;;
    --rt_method) rt_method="$2"; shift 2 ;;
    --RTbound_min_ms) RTbound_min_ms="$2"; shift 2 ;;
    --RTbound_max_ms) RTbound_max_ms="$2"; shift 2 ;;
    --exclude_file) exclude_file="$2"; shift 2 ;;
    --ic_method) ic_method="$2"; shift 2 ;;
    --dry-run) dry_run=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

# Check required arguments
if [[ -z "${model-}" || -z "${cohort-}" ]]; then
  echo "Error: --model and --cohort are required."
  usage
  exit 1
fi

# Build command to run R script
cmd=("Rscript" "./scripts/ppc/run_ppc_pipeline.R") # adjust path if needed

cmd+=("--model" "$model")
cmd+=("--cohort" "$cohort")
cmd+=("--task" "$task")

[[ -n "$ses" ]] && cmd+=("--ses" "$ses")
cmd+=("--group" "$group")
cmd+=("--group_name" "$group_name")
[[ -n "$fit_file" ]] && cmd+=("--fit_file" "$fit_file")
cmd+=("--n_sims" "$n_sims")
cmd+=("--block_size" "$block_size")
[[ -n "$output_dir" ]] && cmd+=("--output_dir" "$output_dir")

if $parallel; then
  cmd+=("--parallel" "--n_cores" "$n_cores")
fi

cmd+=("--steps" "$steps")

if $force; then
  cmd+=("--force")
fi

if $render; then
  cmd+=("--render")
fi

cmd+=("--sampling" "$sampling")
cmd+=("--width_control" "$width_control")
cmd+=("--rt_method" "$rt_method")
cmd+=("--RTbound_min_ms" "$RTbound_min_ms")
cmd+=("--RTbound_max_ms" "$RTbound_max_ms")

[[ -n "$exclude_file" ]] && cmd+=("--exclude_file" "$exclude_file")
cmd+=("--ic_method" "$ic_method")

# Dry run mode
if $dry_run; then
  echo "Dry run: the command would be:"
  printf '%q ' "${cmd[@]}"
  echo
  exit 0
fi

# Run the command
"${cmd[@]}"
