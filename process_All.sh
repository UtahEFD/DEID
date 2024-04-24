#!/bin/bash

# Path to the MATLAB script you want to run
matlab_script="/uufs/chpc.utah.edu/common/home/snowflake3/DEID/Data_Processing_Code/DEID_Processor.m"

# Directory containing folders
base_dir="/uufs/chpc.utah.edu/common/home/snowflake3/DEID/Atwater/JAN"

# Iterate over each folder in the directory
for folder in "$base_dir"/*; do
    # Check if the item is a directory
    if [ -d "$folder" ]; then
        echo "Processing folder: $folder"

        # Navigate to the folder
        cd "$folder" || { echo "Cannot cd into $folder"; exit 1; }

        # Run MATLAB script
        matlab -nodisplay -r "try, run('$matlab_script'), catch ME, disp(getReport(ME)), end, exit"

        # Navigate back to the base directory
        cd "$base_dir" || { echo "Cannot cd back to $base_dir"; exit 1; }
    fi
done

echo "Finished processing all folders."
