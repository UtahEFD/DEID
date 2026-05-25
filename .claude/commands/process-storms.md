# Process All Storms

Process all 16 storms sequentially by editing `working_dir` and `output_dir` directly in `/main/run_deid_processing.m`, running MATLAB, and verifying success before moving to the next storm.

## Storm List

| storm_name | working_dir | output_dir |
|---|---|---|
| mar0723_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/CLN/mar07_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/mar0723_storm |
| mar2023_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/CLN/mar20_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/mar2023_storm |
| mar2423_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/CLN/mar24_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/mar2423_storm |
| dec0423_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/DEC/dec04_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/dec0423_storm |
| dec0923_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/DEC/dec09_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/dec0923_storm |
| jan0524_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/JAN/jan05_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/jan0524_storm |
| jan0824_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/JAN/jan08_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/jan0824_storm |
| jan1024_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/JAN/jan10_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/jan1024_storm |
| feb0624_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/FEB/feb06_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/feb0624_storm |
| feb0724_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/FEB/feb07_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/feb0724_storm |
| feb0824_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/FEB/feb08_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/feb0824_storm |
| feb2024_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/FEB/feb20_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/feb2024_storm |
| dec3024_storm | /uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025/DEC/dec30_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/dec3024_storm |
| jan0225_storm | /uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025/JAN/jan02_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/jan0225_storm |
| feb1425_storm | /uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025/FEB/feb14_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/feb1425_storm |
| apr0125_storm | /uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025/MAR/apr01_storm | /uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/apr0125_storm |

## Steps

For each storm in the list above, perform the following steps in order. Do not move to the next storm until the current one has been verified as successful.

### 1. Edit the script

In `/main/run_deid_processing.m`, find the lines that assign `working_dir` and `output_dir` and replace their values with the current storm's paths.

### 2. Run MATLAB

Run MATLAB with no timeout — storms may take longer than 10 minutes to process. Do not impose any timeout limit. Wait as long as needed for MATLAB to finish.

```bash
module load matlab && matlab -batch "run_deid_processing" > /tmp/storm_log.txt 2>&1
```

Do not create temporary copies of the script (e.g. `run_deid_processing_<storm_name>.m`). Always run the original script directly after editing it in place.

### 3. Verify success

Check that all 4 of the following output files exist in the current storm's `output_dir`:
- A file matching: `DEID_unfilteredParticle_*.csv`
- A file matching: `DEID_filteredParticle_*.csv`
- A file matching: `DEID_aviTotals_*.csv`
- A file matching: `DEID_TS_10min_*.csv`

If any file is missing, print the contents of `/tmp/storm_log.txt` and stop immediately.

### 4. Report and continue

Confirm to the user that the storm succeeded (e.g. "✓ mar0723_storm complete — 1 of 16") then move to the next storm.

### 5. Final report

Once all 16 storms are complete, print a summary listing all storms as confirmed successful.
