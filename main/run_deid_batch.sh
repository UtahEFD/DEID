#!/bin/bash
# Batch-process all day subdirectories within each month directory in SEASON_DIR.

SEASON_DIR=/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2023_2024
SCRIPT=/uufs/chpc.utah.edu/common/home/u6022893/Documents/DEID/dataProcessingCode/main/run_deid_processing.m
LOG_DIR=/uufs/chpc.utah.edu/common/home/u6022893/Documents/DEID/batch_runs

YEAR=$(basename "$SEASON_DIR")
LOG_FILE="$LOG_DIR/deid_${YEAR}_all_months_batch.log"
exec > >(tee -a "$LOG_FILE") 2>&1

for month_dir in "$SEASON_DIR"/*/; do
    MONTH=$(basename "$month_dir")
    echo "###################################################"
    echo "Processing month: $MONTH"
    echo "Time: $(date)"
    echo "###################################################"

    for subdir in "$month_dir"/*/; do
        DAY=$(basename "$subdir")
        export DEID_WORKING_DIR="$subdir"

        # Skip directories with no .avi files
        shopt -s nullglob
        avi_files=("$DEID_WORKING_DIR"/*.avi "$DEID_WORKING_DIR"/*.AVI)
        shopt -u nullglob
        if [ ${#avi_files[@]} -eq 0 ]; then
            echo "Skipping $MONTH/$DAY — no .avi files found."
            continue
        fi

        echo "========================================"
        echo "Starting: $MONTH/$DAY"
        echo "Time: $(date)"
        echo "========================================"
        matlab -nodisplay -nosplash -batch "addpath('/uufs/chpc.utah.edu/common/home/u6022893/Documents/DEID/dataProcessingCode/functions'); addpath('/uufs/chpc.utah.edu/common/home/u6022893/Documents/DEID/dataProcessingCode/main'); run('$SCRIPT')"
        echo "Finished: $MONTH/$DAY at $(date)"
    done

    echo "Done with month: $MONTH at $(date)"
done

echo "All months processed."
