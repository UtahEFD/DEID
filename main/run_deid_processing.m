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

working_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/FEB/allFEB';
output_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/processedData/test2';

%% global variables and physical constants

% specifies resampling period:

time_interval = 600;  % seconds
time_step = seconds(time_interval); % datetime step 

% unit conversions:

% 2025 - present:

% pix_to_m_conversion = .01/40; % m per pix 

% 2023 - 2024:

pix_to_m_conversion = 3.1750e-04; % m per pix

pix_to_m2_conversion = pix_to_m_conversion^2; % m^2 per pix^2
c1 = 10^3; % converts meters to millimeters
max_temp = 145; % max temperature set in colorbar on the physical screen of the tir software
max_int = 255; % maximum intensity that the maximum temperature (145) maps to 
int_to_temp_conversion = max_temp/max_int; % we know that 145 C maps to 255 intensity value in greyscale 

% physical constants:

rho_water = 1000; % density of water [kg/m^3]
mu = 1.5*10^-5;   % viscosity of air [kg/m*s] 

% thresholds & filters:

min_thres = 70; % minimum threshold number in image accepted rbg ([0 255]) scale
minimum_hydro_area = 2; % the minimum number of pixels a hydrometeor must contain to be analyzed 
sort_threshold = 20; % this is the RMS threshold between succesive images of snowflakes used in the sortPostitions_v2. dhiraj calibrated this in the lab. 
minimum_drop_life = 0; % minimum number of frames a drop has to be visable to be processed
areaTol = 0; 
SWEfactor_threshold = 1.85; % maximum value of tolerable SWE factor
evapTime_min = 1/15; % minimum time a snowflake has to appear on hotplate to be processed
evapTime_max = 30; % maximum time a snowflake can appear on hotplate to be processed 
noiseThresh = 999; % # of times a centroid can appear in an .avi file to be considered real 

% DEID specific parameters:

colorbar_image_indexes = [1 1 384 288]; % location of colorbar in pixel locations
crop_index = 43; % use this to specify indices to crop out kapton tape
colorbar_kapton_image_indexes = [1 (colorbar_image_indexes(2)+crop_index) 383 (colorbar_image_indexes(4)-crop_index)]; % location of Kapton tape in pixel locations
k_dLv = 0.002; % calibration constant; in paper, thermal conductivity (k) of water. See sect 4.1 in Dihiraj's paper -> (k/d(_eff))/Latent heat of vaporazation [units?]
l_constant = 2.594e06; % latent heat of vaporazation of water, should be a function of tempertaure (Look at Stull textbook) [J/kg]
% Eqn. (13) in Dhiraj's density paper: c = (L_vv) / (L_ff*C_melt)
% c = hf_rho_coeff
hf_rho_coeff = 1.01e05; % [K*s*m^-1]

%% move to working directory and identify video files

cd(working_dir);

% CALL get_sorted_videos: this returns all .avi file in the directory and
% sorts in time order 

[file_names, vid_date, storm_output] = get_sorted_videos(working_dir); 

%% begin DEID video processing

pbp_table_cell = cell(length(file_names),1);
pbp_table_filtered_cell = cell(length(file_names),1);
avi_summary_table_cell = cell(length(file_names),1);

parfor file_i = 1:length(file_names)
    filename = file_names{file_i};
    disp(['Processing File: ', filename])
    
    % CALL process_one_video - the per-video worker function that will call
    % help functions

    [pbp_table, pbp_table_filtered, avi_summary_table] = process_one_video( ...
        filename, working_dir, ...
        colorbar_image_indexes, colorbar_kapton_image_indexes, ...
        min_thres, minimum_hydro_area, sort_threshold, areaTol, ...
        SWEfactor_threshold, evapTime_min, evapTime_max, ...
        pix_to_m_conversion, pix_to_m2_conversion, int_to_temp_conversion, ...
        k_dLv, l_constant, hf_rho_coeff, rho_water, c1);

    pbp_table_cell{file_i} = pbp_table;
    pbp_table_filtered_cell{file_i} = pbp_table_filtered;
    avi_summary_table_cell{file_i} = avi_summary_table;
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
    pbp_table_retimed = retime_pbp_filtered(pbp_table_filtered, time_step, rho_water, hp_area);
else
    pbp_table_retimed = timetable();
end

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
