[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# DEID Video Processing Pipeline

This repository contains a MATLAB-based pipeline for processing `.avi` video files recorded by the Differential Emissivity Imaging Disdrometer (DEID) of hydrometeors and computing Snow Water Equivalent (SWE) using two approaches:

- Frame-by-frame (FBF) method  
- Particle-by-particle (PBP) method  

The code is designed to process multiple video files efficiently using parallel computing (`parfor`).

---

## 📁 Repository Structure

```
repo/
├── main/
│   ├── run_deid_processing.m                 % main script (call function)
│   ├── DEID_Calibrator.m                     % used for calibrating k/d coefficient 
│   └── DEID_AutoTransfer.m                   % performs all the same functions as 'run_deid_processing.m' but for single .avi files recorded in real time
├── functions/
│   ├── append_gap_row_and_summary.m
│   ├── build_avi_summary_table.m
│   ├── fbf_method.m
│   ├── get_deid_params.m                     % call parameters specific to DEID
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
│   ├── timeAveraged_particle_table.csv        % time averaged particle data file; filtered data 
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
output_dir  = 'example_output';

Then run:

run('main/run_deid_processing.m')

This will generate example output files in:

example_output/ 

### 🏃 Running with Your Own Data

Open: 

main/run_deid_processing.m

Set your directories:

working_dir = 'path/to/avi/files';
output_dir  = 'path/to/save/results';

Begin processing code: 

run('main/run_deid_processing.m')

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
    - Time-averaged SWE results

## 📊 Outputs

The pipeline generates:

    - Particle-by-particle tables: 'DEID_unfilteredParticle_YYYY-MM-DD_HH-MM-SS'
    - Filtered particle-by-particle tables: 'DEID_filteredParticle_YYYY-MM-DD_HH-MM-SS.csv'
    - Per-video summary tables: 'DEID_aviTotals_DD-MM-YYYY.csv'
    - Time-averaged data: 'DEID_TS_MMmin_YYYY-MM-DD_HH-MM-SS.csv'

All outputs are saved to the specified output_dir.

## 🧪 Notes

    - Input .avi files are not stored in this repository; located on the University of Utah's Center for High Performance Computing (CHPC)
    - Data can be made available upon request 
    - Output files are saved externally and are not tracked by Git; also located on CHPC and can be made available upon request
    - The legacy/ folder contains older versions for reference

## 👤 Authors

Ben Silberman 
Dhiraj Singh
Travis Morrison
Alex Blackmer
