# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this codebase does

Processes thermal infrared `.avi` videos from a DEID (hot plate) snow instrument. Each video captures individual hydrometeors (snowflakes, rain drops) landing on a heated Kapton-tape plate. The code extracts per-particle microphysical properties and computes Snow Water Equivalent (SWE) using two independent methods.

## Running the pipeline

This is pure MATLAB — there is no build system or package manager. To process a storm:

1. Open `main/run_deid_processing.m` in MATLAB
2. Set `working_dir` (path to `.avi` files) and `output_dir` (where CSVs are saved)
3. Run the script

Requires MATLAB with:
- **Parallel Computing Toolbox** (`parfor` loop over videos)
- **Image Processing Toolbox** (`regionprops`, `imcrop`, `bwareaopen`, `imfill`)

There is no test suite.

## Architecture

### Entry point
`main/run_deid_processing.m` — loads parameter structs, calls `get_sorted_videos`, runs a `parfor` loop over all `.avi` files calling `process_one_video`, then concatenates results and writes four output CSVs.

### Parameter configuration (edit these for new deployments)
- `functions/get_physical_constants.m` — pixel-to-meter scale (`mPerPix`), temperature/intensity conversion, physical constants (densities, etc.). **`mPerPix` has date-range comments and must be updated when the instrument optics change.**
- `functions/get_thresholds.m` — image segmentation thresholds (`min_thres`, `minimum_hydro_area`) and particle filter bounds (`evapTime_min/max`, `SWEfactor_threshold`).
- `functions/get_deid_params.m` — instrument crop indices (`kapton_indexes`, `colorbar_kapton_image_indexes`) and calibration constants (`k_dLv`, `hf_rho_coeff` from Singh et al. 2024).

### Per-video processing: `process_one_video.m`
Orchestrates four sequential steps for one `.avi`:

1. **`fbf_method.m`** (Frame-By-Frame) — iterates every frame, thresholds the grayscale image to isolate hydrometeors via `regionprops`, computes `∑(area × ΔT)` per frame, integrates to get bulk SWE and a `h_data_cells` array of per-frame particle properties.

2. **`sort_h_data_cells.m` → `sortPositions_v2.m`** — reorganizes the per-frame particle data into tracked particle trajectories (rows = particle IDs, columns = frames).

3. **`pbp_method.m`** (Particle-By-Particle) — iterates over tracked particles, extracts each particle's lifetime properties (mass, area, perimeter, density via heat-flux method, SDI, complexity, shear strength), applies a SWE scaling factor against the FBF result, and returns a timetable where each row is one particle event.

4. **`append_gap_row_and_summary.m` → `build_avi_summary_table.m`** — filters the PBP table (removes residue flags and out-of-range evaporation times), appends a gap-fill summary row at the end of each AVI, and builds a one-row AVI summary timetable.

### Post-processing (back in `run_deid_processing.m`)
- `vertcat` all per-video timetables, `sortrows` by time, recompute cumulative accumulations
- `retime_pbp_filtered.m` — bins filtered PBP data to a regular time grid (default 10 min), computing mean/sum columns, density, shear strength, SWE, and snow depth per bin

### Output files (written to `output_dir`)
| File | Contents |
|------|----------|
| `DEID_unfilteredParticle_*.csv` | All detected particles (timetable, one row per particle) |
| `DEID_filteredParticle_*.csv` | Filtered particles (residue and evap-time filters applied) |
| `DEID_aviTotals_*.csv` | One summary row per `.avi` file |
| `DEID_TS_10min_*.csv` | 10-minute time-averaged timetable |

## Key domain concepts

- **FBF SWE**: bulk SWE from integrating thermal signal across all pixels in each frame — more stable but no per-particle info.
- **PBP SWE**: SWE summed from individual particle masses — then scaled by `SWE_factor` (ratio of FBF to PBP totals, capped at `SWEfactor_threshold`) to reconcile the two methods.
- **SDI** (Snowflake Disequilibrium Index): ratio of actual snowflake area to equivalent water droplet area — measures crystal complexity.
- **Complexity (Cx)**: bounding-rectangle area / actual area — another shape metric.
- **Shear Strength**: `σ_ice × SDI × (ρ/ρ_ice)^Cx` — used as a snow strength proxy.
- **Kapton tape**: the heating element strip whose maximum pixel intensity gives the plate temperature each frame.
- **`h_mass_fbf_min`**: minimum per-frame FBF mass (baseline noise floor), subtracted before SWE integration.

## Legacy code
`legacy/` contains older standalone scripts (`DEID_Processor.m`, Dhiraj's scripts) kept for reference. Do not modify them; all active development is in `functions/` and `main/`.
