[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
## 📄 License

This project is licensed under the MIT License – see the LICENSE file for details.
# DEID Video Processing Pipeline

This repository contains a MATLAB-based pipeline for processing `.avi` video files recorded by the Differential Emissivity Imaging Disdrometer (DEID). Using the DEID, we measure hydrometeor properties including density, area, and diameter, while collecting storm properties such as Snow Water Equivalence (SWE) and Precipitation Intensity (PI). We use two approaches:

- Frame-by-frame (FBF) method
- Particle-by-particle (PBP) method

The FBF method analyzes `.avi` files frame-by-frame. For each frame, all hydrometeors are identified and quantified collectively. The PBP method tracks and quantifies individual hydrometeors across all frames where they are present. The FBF method provides the most accurate SWE measurements, while the PBP method provides accurate particle density measurements.

The code is designed to process multiple video files efficiently using parallel computing (`parfor`).

---

## 📁 Repository Structure

```
repo/
├── main/
│   ├── run_deid_processing.m                 % main script (call function)
│   ├── run_deid_batch.sh                     % bash script to batch-process multiple subdirectories
│   └── DEID_Calibrator.m                     % used for calibrating k/d coefficient
├── functions/
│   ├── append_gap_row_and_summary.m
│   ├── build_avi_summary_table.m
│   ├── fbf_method.m
│   ├── get_deid_params.m                     % load parameters specific to DEID
│   ├── get_physical_constants.m              % call physical constants and any unit conversions
│   ├── get_sorted_videos.m
│   ├── get_thresholds.m                      % call thresholds used for filtering and cleaning data
│   ├── pbp_method.m
│   ├── process_one_video.m
│   ├── retime_pbp_filtered.m
│   ├── sort_h_data_cells.m
│   └── sortPositions_v2.m
├── legacy/
│   ├── old_script/                            % previous versions of scripts
│   └── dhiraj_script/                         % original code developed by Dhiraj Singh
├── example_data/
│   └── DEID_sampleVideo.avi                   % a cropped (~35MB) video file used for testing
├── example_output/
│   ├── unfiltered_particle_table.csv          % unfiltered particle-by-particle data file
│   ├── filtered_particle_table.csv            % filtered particle-by-particle data file
│   ├── timeAveraged_particle_table.csv        % time-averaged particle data file; filtered data
│   └── summary_table.csv                      % appended summary of each avi file
├── README.md
└── .gitignore
```
---

## 🚀 Getting Started

### Requirements

- MATLAB (R2021a or newer recommended)
- Image Processing Toolbox
- Parallel Computing Toolbox (for `parfor`)

---

### ⚡ Quick Start (Using Sample Data)

A small example video is included in the repository:

example_data/DEID_sampleVideo.avi

Set your directories:

working_dir = 'example_data';
output_dir  = working_dir;

Then run:

run('main/run_deid_processing.m')

Outputs will be saved to the same directory as the input data.

### 🏃 Running with Your Own Data

#### Single directory (MATLAB)

Open `main/run_deid_processing.m` and set `working_dir` near the top:

```matlab
working_dir = 'path/to/avi/files';
```

Output files are saved to the same directory as the input (`output_dir = working_dir`).

The script auto-detects the sensor from the path: the path must contain either `"snowflake3"` or `"snowflake4"`. You can also override `working_dir` without editing the script by setting the environment variable `DEID_WORKING_DIR` before launching MATLAB.

Begin processing:

```matlab
run('main/run_deid_processing.m')
```

#### Batch processing multiple subdirectories (Bash)

Use `main/run_deid_batch.sh` to process all subdirectories within a base directory sequentially. Edit the top of the script to set:

```bash
BASE=/path/to/parent/directory    # all immediate subdirectories will be processed
LOG_DIR=/path/to/log/directory
```

The script iterates over every subdirectory in `BASE`, sets `DEID_WORKING_DIR` for each, invokes MATLAB, and writes a log to `LOG_DIR/deid_<month>_batch.log`.

Run with:

```bash
bash main/run_deid_batch.sh
```

Or in the background:

```bash
nohup bash main/run_deid_batch.sh &
```

## ⚙️ What the Code Does

### For each .avi file, the pipeline:

1. Frame-by-frame processing
    - Converts frames to grayscale
    - Identifies hydrometeors
    - Computes area–temperature products
    - Calculates SWE (FBF method)
2. Tracking and sorting
    - Matches hydrometeors across frames
    - Organizes data into consistent structures
3. Particle-by-particle analysis
    - Tracks individual hydrometeors through time
    - Computes mass, density, evaporation time, and SWE contribution
4. Filtering and corrections
    - Removes noisy or non-physical particles
    - Applies SWE correction factor
5. Output generation
    - Particle-level data tables
    - Filtered datasets
    - Per-video summary tables
    - Time-averaged SWE results (10-minute bins)

## 📊 Outputs

All outputs are saved to the same directory as the input `.avi` files (`output_dir = working_dir`).

The pipeline generates four files per storm:

| File | Contents |
|------|----------|
| `DEID_unfilteredParticle_YYYY-MM-DD_HH-MM-SS.csv` | All detected particles (one row per particle) |
| `DEID_filteredParticle_YYYY-MM-DD_HH-MM-SS.csv` | Filtered particles (residue and evap-time filters applied) |
| `DEID_aviTotals_<storm>.csv` | One summary row per `.avi` file |
| `DEID_TS_10min_YYYY-MM-DD_HH-MM-SS.csv` | 10-minute time-averaged timetable |

## 🧪 Notes

    - Input `.avi` files are stored on the University of Utah CHPC and are not included in this repository
    - Data and outputs can be made available upon request
    - Output files are not tracked by Git
    - The `legacy/` folder contains older versions for reference

## 👤 Authors

Ben Silberman<br>
Dhiraj Singh<br>
Travis Morrison<br>
Alex Blackmer<br>
