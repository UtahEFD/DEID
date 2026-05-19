%% DEID AVI File Processing Code
%
% OUTPUTS (for each .avi file processsed):
%
%   - filtered particle-by-particle (pbp) .csv file: DEID_filteredParticle_YYYY-MM-DD_HH-mm-ss.csv
%   - unfiltered particle-by-particle (pbp) .csv file: DEID_unfilteredParticle_YYYY-MM-DD_HH-mm-ss.csv
%   - an appended .csv file containing a row of data summarizing each .avi
%       file processed; 'DEID_aviTotals.csv'
%
% data is processed using a parallelized forloop (PARFOR)

% AUTHOR : Dhiraj Singh, Benjamin Silberman, Travis Morrison, Alex Blackmer
                                       
clear, clc, close all

%% set filepath, output directory, and file name for saving  

working_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/CLN/mar07_storm';
output_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/stormData/mar0723_storm';

%% load parameter structures

phys   = get_physical_constants();
thresh = get_thresholds();
deid   = get_deid_params();

%% move to working directory and identify video files

cd(working_dir);

% CALL get_sorted_videos: this returns all .avi file in the directory and
% sorts in time order 

[file_names, vid_date, storm_output] = get_sorted_videos(working_dir); 

%% begin DEID video processing

% specify resampling period for time averaged data:

time_step = seconds(600);

% initialize cells for parfor loop: 

pbp_table_cell = cell(length(file_names),1);
pbp_table_filtered_cell = cell(length(file_names),1);
avi_summary_table_cell = cell(length(file_names),1);
plate_temp_cell = cell(length(file_names),1);
back_temp_cell = cell(length(file_names),1);
k_dLv_cell = cell(length(file_names),1);
time_fbf_cell = cell(length(file_names),1);

parfor file_i = 1:length(file_names)
    filename = file_names{file_i};
    disp(['Processing File: ', filename])

    % CALL process_one_video - the per-video worker function that will call
    % help functions

    [pbp_table, pbp_table_filtered, avi_summary_table, plate_temp_fbf, back_temp_fbf, k_dLv_calibrate_fbf, time_series_fbf] = process_one_video( ...
        filename, working_dir, deid.kapton_indexes, ...
        deid.colorbar_image_indexes, deid.colorbar_kapton_image_indexes, ...
        thresh.min_thres, thresh.minimum_hydro_area, thresh.sort_threshold, thresh.areaTol, ...
        thresh.SWEfactor_threshold, thresh.evapTime_min, thresh.evapTime_max, ...
        phys.mPerPix, phys.m2PerPix2, phys.int_to_temp_conversion, ...
        deid.k_dLv, deid.hf_rho_coeff, phys.rho_water, phys.rho_ice, phys.sigma_ice, ...
        phys.mmPerM);

    pbp_table_cell{file_i} = pbp_table;
    pbp_table_filtered_cell{file_i} = pbp_table_filtered;
    avi_summary_table_cell{file_i} = avi_summary_table;
    plate_temp_cell{file_i} = plate_temp_fbf;
    back_temp_cell{file_i} = back_temp_fbf;
    k_dLv_cell{file_i} = k_dLv_calibrate_fbf;
    time_fbf_cell{file_i} = time_series_fbf;
end

%% now put back into one large time table:

pbp_table = vertcat(pbp_table_cell{:});
pbp_table = sortrows(pbp_table, 'Time');
pbp_table.("PBP SWE Accumulation (mm)") = cumsum(pbp_table.("PBP SWE (mm)"));
pbp_table.("FBF SWE Accumulation (mm)") = cumsum(pbp_table.("FBF SWE (mm)"));
pbp_table.("PBP Snow Accumulation (mm)") = cumsum(pbp_table.("PBP Snow (mm)"));
pbp_table.("FBF Snow Accumulation (mm)") = cumsum(pbp_table.("FBF Snow (mm)"));

pbp_table_filtered = vertcat(pbp_table_filtered_cell{:});
pbp_table_filtered = sortrows(pbp_table_filtered, 'Time');
pbp_table_filtered.("PBP SWE Accumulation (mm)") = cumsum(pbp_table_filtered.("PBP SWE (mm)"));
pbp_table_filtered.("FBF SWE Accumulation (mm)") = cumsum(pbp_table_filtered.("FBF SWE (mm)"));
pbp_table_filtered.("PBP Snow Accumulation (mm)") = cumsum(pbp_table_filtered.("PBP Snow (mm)"));
pbp_table_filtered.("FBF Snow Accumulation (mm)") = cumsum(pbp_table_filtered.("FBF Snow (mm)"));

avi_summary_table = vertcat(avi_summary_table_cell{:});
avi_summary_table = sortrows(avi_summary_table, 'Time');
hp_area = avi_summary_table.("Hot Plate Area")(1);

%% time average data here:

if ~isempty(pbp_table_filtered)
    pbp_table_retimed = retime_pbp_filtered(pbp_table_filtered, time_step, phys.rho_water, phys.rho_ice, phys.sigma_ice, hp_area);
else
    pbp_table_retimed = timetable();
end

%% storm-level FBF diagnostic plot:

all_time_fbf   = [time_fbf_cell{:}];
all_plate_temp = vertcat(plate_temp_cell{:});
all_back_temp  = vertcat(back_temp_cell{:});
all_k_dLv      = vertcat(k_dLv_cell{:});

[all_time_fbf, sort_idx] = sort(all_time_fbf);
all_plate_temp = all_plate_temp(sort_idx);
all_back_temp  = all_back_temp(sort_idx);
all_k_dLv      = all_k_dLv(sort_idx);

fig = figure('Visible', 'off');

subplot(3,1,1)
plot(all_time_fbf, all_back_temp, 'b.', 'MarkerSize', 3)
ylabel('Temp (°C)')
title('Background Temperature')
grid on

subplot(3,1,2)
plot(all_time_fbf, all_plate_temp, 'r.', 'MarkerSize', 3)
ylabel('Temp (°C)')
title('Plate Temperature')
grid on

subplot(3,1,3)
plot(all_time_fbf, all_k_dLv, 'k.', 'MarkerSize', 3)
ylabel('k/d Coefficient')
xlabel('Time')
title('K/D Coefficient')
grid on

sgtitle(['Per-Frame Summary: ', storm_output], 'Interpreter', 'none')
saveas(fig, fullfile(output_dir, ['DEID_fbf_diagnostics_', storm_output, '.png']));
close(fig);

%% save processed tables:

if ~isempty(pbp_table_retimed)
    startTime = datestr(pbp_table_retimed.Time(1), 'yyyy-mm-dd_HH-MM-ss');
else
    startTime = datestr(now, 'yyyy-mm-dd_HH-MM-ss');
end

writetimetable(pbp_table, [output_dir, '/DEID_unfilteredParticle_', startTime, '.csv']);
writetimetable(pbp_table_filtered, [output_dir, '/DEID_filteredParticle_', startTime, '.csv']);
writetimetable(avi_summary_table, [output_dir, '/DEID_aviTotals_', storm_output, '.csv']);
if ~isempty(pbp_table_retimed)
    writetimetable(pbp_table_retimed, [output_dir,'/DEID_TS_10min_', startTime, '.csv']);
end

[~, parent_dir, ~] = fileparts(pwd);
disp(['Saved Output for: ', parent_dir])
