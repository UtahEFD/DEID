#!/bin/bash
# Retry script: process only specific January subdirectories that were missed/failed.

BASE=/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/2023_2024/Atwater/JAN
SCRIPT=/uufs/chpc.utah.edu/common/home/u6022893/Documents/DEID/dataProcessingCode/main/run_deid_processing.m
LOG_DIR=/uufs/chpc.utah.edu/common/home/u6022893/Documents/DEID/batch_runs

MONTH=$(basename "$BASE" | tr '[:upper:]' '[:lower:]')
YEAR=$(basename "$(dirname "$(dirname "$BASE")")")
LOG_FILE="$LOG_DIR/deid_${YEAR}_${MONTH}_retry.log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Only process these specific subdirectories
TARGETS=("jan15" "jan24" "jan25" "jan26" "jan27")

for subdir in "${TARGETS[@]}"; do
    export DEID_WORKING_DIR="$BASE/$subdir"

    if [ ! -d "$DEID_WORKING_DIR" ]; then
        echo "Skipping $subdir — directory not found."
        continue
    fi

    # Skip directories with no .avi files
    shopt -s nullglob
    avi_files=("$DEID_WORKING_DIR"/*.avi "$DEID_WORKING_DIR"/*.AVI)
    shopt -u nullglob
    if [ ${#avi_files[@]} -eq 0 ]; then
        echo "Skipping $subdir — no .avi files found."
        continue
    fi

    echo "========================================"
    echo "Starting: $DEID_WORKING_DIR"
    echo "Time: $(date)"
    echo "========================================"
    matlab -nodisplay -nosplash -batch "addpath('/uufs/chpc.utah.edu/common/home/u6022893/Documents/DEID/dataProcessingCode/functions'); addpath('/uufs/chpc.utah.edu/common/home/u6022893/Documents/DEID/dataProcessingCode/main'); run('$SCRIPT')"
    echo "Finished: $subdir at $(date)"
done

echo "All targeted subdirectories processed."
