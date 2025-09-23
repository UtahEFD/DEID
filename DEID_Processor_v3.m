%% DEID AVI File Processing Code
% Outputs particle by particle csclear,v file to be post processed by
% DEID_AVI_Processor.m
% AUTHOR : Dhiraj Singh, Benjamin Silberman, Travis Morrison, Alex Blackmer
                                       
clear, clc, close all
%% set filepath, output directory, and file name for saving  
working_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/DEC/DEC_all';
output_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/test';
storm_output = '_dec042023_storm';

%% global variables and physical constants
% specifies resampling period:
time_interval = 600;  % seconds
time_step = seconds(time_interval); % datetime step 
% unit conversions:
mm_to_inches = 1/25.4; % [mm/in]
pix_to_m_conversion = 3.1750e-04; % m per pix 
pix_to_m2_conversion = pix_to_m_conversion^2; % m^2 per pix^2
c1 = 10^3; % converts meters to millimeters
c2 = 3.6e06; % converts m/s to mm/hr [(mm/hr)*(s/m)]
hour_in_seconds = 3600;
% physical constants:
rho_water = 1000; % density of water [kg/m^3]
mu = 1.5*10^-5;   % viscosity of air [kg/m*s] 
% thresholds & filters:
flagTolerance = 1e-6; % how much the area and temp diff of a given snowflake over its life cycle can change for it to be treated as resiude
SWEfactor_threshold = 1.5; % maximum value of tolerable SWE factor
min_thres = 70; % minimum threshold number in image accepted rbg ([0 255]) scale
residue_filter = 0.005; % max weight of snowflake to process [kg]
evapTime_min = 1/15; % minimum time a snowflake has to appear on hotplate to be processed
evapTime_max = 30; % maximum time a snowflake can appear on hotplate to be processed
minimum_drop_life = 0; % minimum number of frames a drop has to be visable to be processed
colorbar_max_temp = 145; % max temperature set in colorbar on the physical screen of the tir software
sort_threshold = 20; % this is the RMS threshold between succesive images of snowflakes used in the sortPostitions_v2. dhiraj calibrated this in the lab. 
% DEID specific parameters:
colorbar_image_indexes = [1 1 384 288]; % location of colorbar in pixel locations
crop_index = 43; % use this to specify indices to crop out kapton tape
colorbar_kapton_image_indexes = [1 (colorbar_image_indexes(2)+crop_index) 383 (colorbar_image_indexes(4)-crop_index)]; % location of Kapton tape in pixel locations
k_dLv = 0.0030; % calibration constant; in paper, thermal conductivity (k) of water. See sect 4.1 in Dihiraj's paper -> (k/d(_eff))/Latent heat of vaporazation [units?]
l_constant = 2.594e06; % latent heat of vaporazation of water, should be a function of tempertaure (Look at Stull textbook) [J/kg]
% Eqn. (13) in Dhiraj's density paper: c = (L_vv) / (L_ff*C_melt)
% c = hf_rho_coeff
hf_rho_coeff = 1.01e05; % [K*s*m^-1]
%% move to working directory and identify video files 
% sets path for .m to run in:
cd(working_dir) 
% get a list of all files and folders in this directory
directory = dir(".");
% initialize an empty cell array to store the file names
file_names = cell(1, length(directory));
count = 0; 
% loop through each item in the input directory
for file_i = 1:length(directory)
    % check if the item is a file (not a folder) and if it ends with .avi
    if ~directory(file_i).isdir && endsWith(directory(file_i).name, '.avi', 'IgnoreCase', true)
        count = count +1; 
        file_names{count} = directory(file_i).name;
    end
end
file_names = file_names(1:count);
% **when testing**
file_names = file_names(1:4);

%% initialize output tables
% time series output table
ts_col_names = {'Time', 'Complexity', 'SDI', 'Mass', 'Volume', 'Diameter', 'Surface Area', 'Void Space', 'Density [kg/m^3]', 'SWE [mm]','Snow [mm]'};
ts_col_types = {'datetime', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
ts_output_table = table('Size', [0, length(ts_col_names)], ...
                         'VariableNames', ts_col_names, ...
                         'VariableTypes', ts_col_types);
% DEID summary table 
% summary_col_names = {'Time', 'Duration', 'Complexity', 'SDI', 'HFD Density (kg*m^-3)', 'PBP SWE (mm)', 'FBF SWE (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'Hot Plate Area', 'SWE Factor', 'Min FBF Mass'};
% summary_col_types = {'datetime', 'datetime', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
% avi_summary_table = table('Size', [1, length(summary_col_names)], ...
                         % 'VariableNames', summary_col_names, ...
                         % 'VariableTypes', summary_col_types);
 
%% begin DEID video processing

% Parfor loop parallelizes processing by distributing each video file to a
% Matlab worker on each CPU core. 

% preallocate variables:
pbp_table_cell = cell(length(file_names),1);
pbp_table_filtered_cell = cell(length(file_names),1); 
avi_summary_table_cell = cell(length(file_names),1);

 parfor file_i = 1:length(file_names)
    filename = file_names{file_i};
    disp(['Processing File: ', filename])
    vid=VideoReader(filename);
    
    % get metaData for video processing 
    vid_dir = dir(filename);
    vid_fps = vid.FrameRate;
    num_frames = vid.NumFrames;

    % gets video start and end date to construct time series
    vid_end_time = datetime([vid_dir.date]);
    
    % create a time series of date times starting from (video end time - video duration) and ending at the video end time
    time_series = datetime(vid_end_time - (0:num_frames) * seconds(1/vid_fps), 'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
    time_series = flip(time_series);  % Flips time series to be chronologically ordered
    vid_start_time = datetime(time_series(1));

    % preallocate cell for hydrometor data
    h_data = cell(num_frames,1); 

    %% "frame by frame method"; this is how we obtain SWE for each .avi file 
    % preallocate variables saved in loop for speed:
    plate_temp = nan(num_frames,1);
    sum_h_area_times_dt = nan(num_frames,1);
    % initialize videoWriter for this file:
    [~, name, ~] = fileparts(filename); 
    filter_filename = fullfile(output_dir, [name '_filtered.avi']); 
    filter_avi = VideoWriter(filter_filename);
    filter_avi.FrameRate = vid.FrameRate; 
    open(filter_avi); % without opening the file for writing, we can't write any new frames 
    % enter loop to process images: 
    for frame_ii = 1:num_frames     
        frame = read(vid, frame_ii);  
        frame_gray = im2gray(frame); % convert frame of interest to gray scale
        frame_gray_cropped_wKapton = imcrop(frame_gray, colorbar_image_indexes);% crop out colorbar
        plate_temp(frame_ii) = max(max(double(frame_gray_cropped_wKapton))); % this assumes max temperature in image is the plate temperature with Kapton tape 
        frame_gray_cropped = imcrop(frame_gray, colorbar_kapton_image_indexes); % back to orginal grayscale image... now remove colorbar and kapton tape from image
        frame_filtered = frame_gray_cropped > min_thres; % removed below min threshold, on rbg ([0, 255]) scale 
        frame_filtered_filled = imfill(frame_filtered, 'Holes'); % clean up Hydrometeors
        % imshow(frame_filtered_filled)
        % convert frame_filtered_filled (a logical array) to uint8 for
        % grayscale video writing: 
        filter_frameOut = uint8(frame_filtered_filled) * 255;
        writeVideo(filter_avi, filter_frameOut);
        % get hydrometeor properties: 
        h_geo_prop = regionprops(frame_filtered_filled, 'MajorAxisLength', 'MinorAxisLength', 'Centroid', 'Area','BoundingBox'); % returns the centroid, the area of each blob, and the bounding box (left, top, width, height).
        % if no properties are found, go to next frame: 
        if (isempty(h_geo_prop))
            continue;
        end
        % build hydrometeor property matrices from regionprops values: 
        h_bounding_box = cat(1,h_geo_prop.BoundingBox); % concat all values to build matrix of bounding box indices
        h_centroid = round(cat(1, h_geo_prop.Centroid)); % concat all values to build matrix of centroids
        h_major = cat(1,h_geo_prop.MajorAxisLength);
        h_minor = cat(1,h_geo_prop.MinorAxisLength);
        h_area_pix = cat(1, h_geo_prop.Area); % concat all values to build matrix of Hydrometeor areas in pixels
        % convert hydrometeor areas to m^2 and lengths to m: 
        h_area = h_area_pix .* pix_to_m2_conversion; 
        h_major_axis = h_major * pix_to_m_conversion;
        h_minor_axis = h_minor * pix_to_m_conversion;
        % get the hydrometeor ellipse area: 
        h_elipse_area = (pi*h_major_axis .* h_minor_axis)/4;
        % difference in temperature of each centroid and the plate:
        h_centroid_i = sub2ind(size(frame_gray_cropped), h_centroid(:, 2), h_centroid(:, 1)); % find the linear index of the centriods in orginal image
        h_centroid_values = double(frame_gray_cropped(h_centroid_i)); % intensities of centroid pixels of snow
        plate_h_dtemp = colorbar_max_temp - (h_centroid_values .* (colorbar_max_temp / plate_temp(frame_ii)));
        % product of hydrometeor area with the temp difference:
        h_area_times_dtemp = h_area .* plate_h_dtemp;         
        % sum the product of individual area and temp. diff in each frame:
        % **this is how we obtain hydrometeor mass using fbf method**
        sum_h_area_times_dt(frame_ii)=sum(h_area_times_dtemp);         
        % build large matrix of Hydrometeor data
        h_data{frame_ii} = cat(2, h_centroid, h_area, plate_h_dtemp, h_elipse_area, h_major_axis, h_minor_axis); 
    end
    close(filter_avi); % closes and effectivley saves the filtered .avi file
    % frame by frame SWE calculation
    hp_area = size(frame_gray_cropped,1) * size(frame_gray_cropped,2) * pix_to_m2_conversion; % hotplate area     
    h_mass_fbf = (k_dLv*sum_h_area_times_dt) / vid_fps; % total mass evaporates in each frame
    h_mass_fbf_min = min(h_mass_fbf); % we know the plate should be empty when it is not snowing..
    h_mass_fbf = h_mass_fbf - h_mass_fbf_min; % subtract off min mass on a frame to account for any resiude
    SWE_fbf = h_mass_fbf / hp_area; 
    % find the minimum SWE in all frames within a video, and subtract from
    % SWE (way of handling residue) 
    fbf_SWE_min = min(SWE_fbf(SWE_fbf ~=0));
    SWE_fbf = SWE_fbf - fbf_SWE_min; 
    time_series_fbf = time_series(1:length(SWE_fbf));
    
    %% call sortPositions_v2.m to place snowflakes in the same row across 
    % multiple frames, making tracking possible over time
    % create new cell array for sorted data: 
    h_data_sorted = cell(size(h_data));
    h_data_sorted{1} = h_data{1};
    % loop over every data cell corresponding to each frame:
    % Uses "sortPositions_v2.m" - a code that Dhiraj wrote 
    for frame_jj = 2:num_frames
        h_data_sorted{frame_jj} = h_data{frame_jj};
        h_data_sorted{frame_jj} = sortPositions_v2(h_data_sorted{frame_jj-1}, h_data_sorted{frame_jj}, sort_threshold);
    end

    % return frame with max number of Hydrometeors 
    max_h_obs = max(cellfun(@(x) size(x, 1), h_data_sorted, 'UniformOutput', 1));
    % now pad the data with zeros so the frames all have the same number:  
    h_data_sorted = cellfun(@(x) cat(1, x, zeros(max_h_obs - size(x, 1), 7)),...
        h_data_sorted, 'UniformOutput', 0);
    
    %% isolating the variables and put them into a matrix to work with
    
    % For reference: [Hydrometeor_Centroid, Hydrometeor_Area,
    % Plate_Hydrometeor_DeltaT, Hydrometeor_ellipse_area,
    % Hydrometeor_Major_Axis,Hydrometeor_Minor_Axis]
    
    h_area_final = cellfun(@(x) x(:, 3), h_data_sorted, 'UniformOutput', 0);
    dT_final = cellfun(@(x) x(:, 4), h_data_sorted, 'UniformOutput', 0);
    ellipse_area_final = cellfun(@(x) x(:, 5), h_data_sorted, 'UniformOutput', 0);
    maj_axis_final = cellfun(@(x) x(:, 6), h_data_sorted, 'UniformOutput', 0);
    min_axis_final = cellfun(@(x) x(:, 7), h_data_sorted, 'UniformOutput', 0);
    
    % convert to matrix: 
    h_area_final = cat(2,h_area_final{:});
    dT_final = cat(2,dT_final{:});
    ellipse_area_final = cat(2,ellipse_area_final{:});
    maj_axis_final = cat(2,maj_axis_final{:});
    min_axis_final = cat(2,min_axis_final{:});

    %% "particle by particle method" 
    % this method uses the data obtained from fbf, but now isolates each
    % snowflake per frame 
    all_h_appears_ind = []; % initial time index of each snowflake in time array
    h_area = {}; % hydrometeor area over time
    h_delta_temp = {}; % hydrometeor delta temp over time 
    hArea_range = []; % difference between max and min values of hydrometeor area
    deltaTemp_range = []; % difference between max and min values of hydrometeor deltaTemp
    deltaTemp_residue_flags = []; % flags for when deltaTemp does not change
    hArea_residue_flags = []; % flags for when area does not change 
    h_max_area = {}; % hydrometeor maximum area
    h_ellipse_area = {}; % hydrometeor Max circumscribed area
    h_max_maj_axis = {}; % hydrometeor Max major axis
    h_max_min_axis = {}; % hydrometeor Max minor axis
    h_delta_time = {}; % hydrometeor time to melt and evaporate
    h_delta_temp_max = {}; % hydrometeor centroid's max intensity
    h_delta_temp_mean = {}; % hydrometeor centroid's mean intensity
    h_mass_pbp = {}; % hydrometeor calculated mass

    % loop over all Hydrometeors for each frame: 
    % this is where the cleaning occurs!

    for h_ii = 1:max_h_obs
        h_appear_evap_bool = diff(h_area_final(h_ii, :) > 0); %creates an array 
        % of indices that are booleans. +1 is when snowflake appears and -1
        % when it evaps. this is finding the difference in a logical array!

        h_appears_ind = find(h_appear_evap_bool > 0); % finds the index when the hydrometor appears
        if isempty(h_appears_ind)
            continue;
        end
        h_evaps_ind = find(h_appear_evap_bool < 0); % hydrometeor dissapears index
        all_h_appears_ind = cat(2,all_h_appears_ind,h_appears_ind); % concatenate all snowflake arrival indices
        if isempty(h_evaps_ind) % get rid of hydrometeors which do no evaporate completely
            continue; 
        else
        % isolate the hydrometeor properties for its lifetime: 
        h_dT_tmp = dT_final(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
        h_ellipse_area_tmp = ellipse_area_final(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
        h_maj_axis_tmp = maj_axis_final(h_ii,h_appears_ind(1):h_evaps_ind(end)+1);
        h_min_axis_tmp = min_axis_final(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
        % prepare area of each hydrometeor for temporal intergration: 
        h_area_tmp = h_area_final(h_ii, h_appears_ind(1):h_evaps_ind(end)+1); 
        h_area_tmp_bool = h_area_tmp > 0; % isolate for positive values
        h_area_tmp_bool = bwareaopen(h_area_tmp_bool, minimum_drop_life); % any hydrometeor which lives for less than 'minimum_Drop_Life' is discarded 
        propstemp = regionprops(h_area_tmp_bool, 'PixelIdxList'); % finds how many 'continuous' snowflakes there are: indices of contiguous '1s' 
        % this grabs properties of each snowflake over their 'life': 
            for jj = 1:numel(propstemp)
                h_area{end+1} = h_area_tmp(propstemp(jj).PixelIdxList); % snowflake area over time
                h_delta_temp{end+1} = h_dT_tmp(propstemp(jj).PixelIdxList); % snowflake delta temp over time
                h_max_area{end+1} = max(h_area_tmp(propstemp(jj).PixelIdxList)); % max snowflake area over time
                h_ellipse_area{end+1} = max(h_ellipse_area_tmp(propstemp(jj).PixelIdxList)); % max ellipse area over time
                h_max_maj_axis{end+1} = max(h_maj_axis_tmp(propstemp(jj).PixelIdxList)); % max major axis over time
                h_max_min_axis{end+1} = max(h_min_axis_tmp(propstemp(jj).PixelIdxList)); % max minor axis over time 
                h_delta_time{end+1} = numel(propstemp(jj).PixelIdxList); % time of snowflake's life 
                h_delta_temp_max{end+1} = max(h_dT_tmp(propstemp(jj).PixelIdxList)); % max temperature difference between snowflake and plate 
                h_delta_temp_mean{end+1} = mean(h_dT_tmp(propstemp(jj).PixelIdxList)); % mean temperature difference between snowflake and plate
                % filtering:
                hArea_range(end+1) = range(h_area{end}); 
                deltaTemp_range(end+1) = range(h_delta_temp{end}); 
                deltaTemp_residue_flags(end+1) = (range(h_delta_temp{end}) < flagTolerance); % flag snowflakes whose delta temp does not change over time
                hArea_residue_flags(end+1) = (range(h_area{end}) < flagTolerance); % flag snowflakes whose area does not change over time
                % multiply the area by the temperature difference for that hydrometeor:
                h_area_times_dT_pbp = h_area_tmp(propstemp(jj).PixelIdxList) .* h_dT_tmp(propstemp(jj).PixelIdxList);
                % now integrate (sum over snowflake's life):
                h_mass_pbp{end+1} = sum(k_dLv * h_area_times_dT_pbp); % total mass that each snoflake contains across its entire life cylce
            end
        end
    end

    %% convert to matrixes
    h_area = cell2mat(h_area);
    h_delta_temp = cell2mat(h_delta_temp); 
    h_mass_pbp=cell2mat(h_mass_pbp); % hydrometeor mass
    h_max_area=cell2mat(h_max_area); % actual area
    h_ellipse_area=cell2mat(h_ellipse_area); % ellipse area
    h_max_maj_axis=cell2mat(h_max_maj_axis); % height
    h_max_min_axis=cell2mat(h_max_min_axis); % width
    h_delta_temp_max=cell2mat(h_delta_temp_max); % temperature diffrence between plate and water droplet using max intensity
    h_delta_temp_mean=cell2mat(h_delta_temp_mean); % temperature diffrence between plate and water droplet using mean intensity 
    % now multiples mass by "dt":
    h_mass_pbp = h_mass_pbp / vid_fps; 
    h_diam = ((4/pi) * h_max_area).^(1/2); % convert area to diameter of hydrometeor
    h_water_eq_diameter = (6 .* h_mass_pbp / (pi * rho_water)).^(1/3); % water equi. diameter. Not sure where this is from
    h_evap_time = cell2mat(h_delta_time) * (1 / vid_fps); % evaporation time
    h_vol_sph = (3/4) * h_max_area.^(3/2); % spherical volume
    h_rho_sph = h_mass_pbp ./ h_vol_sph; % density calculation: spherical assumption
    h_energy_per_time = h_mass_pbp * l_constant ./ (h_max_area .* h_evap_time); % heat flux method: energy per unit area per time
    h_height = h_evap_time .* h_delta_temp_mean; 
    h_rho_hfd = (hf_rho_coeff * h_mass_pbp) ./ (h_max_area .* h_evap_time .* h_delta_temp_mean); 
    h_vol_hfd = h_mass_pbp ./ h_rho_hfd; % volume of each snowflakes using mean heat flux method density
    h_initial_time_indices = all_h_appears_ind(1:length(h_mass_pbp));
    h_initial_time = time_series(h_initial_time_indices);

    %% conversions and calculations using PBP data
    terminal_vel = terminalVelocity(h_max_area, h_ellipse_area, h_mass_pbp); % terminal velocity 
    void_space = (h_max_area - h_ellipse_area).^(1/2); % void space      
    for i = 1:length(void_space)
        if void_space(i) > 0
            void_space(i) = void_space(i);
        else
            void_space(i) = 0;
        end
    end
    surface_area_eq = ((9*pi)/16)^(1/3) * (h_mass_pbp ./ rho_water).^(2/3); % surface area of equivalent water droplet (See POF Singh et al. 2023)
    complexity = h_ellipse_area ./ h_max_area; % complexity (See CRST Morrison et al. 2023)
    sdi = h_max_area ./ surface_area_eq ; % SDI (See CRST Morrison et al. 2023)
    SWE_pbp = c1 * h_mass_pbp ./ (rho_water * hp_area); % [mm]
    SWE_factor = sum(SWE_fbf) / sum(SWE_pbp); % swe factor is used to adjust pbp SWE 
    if SWE_factor > SWEfactor_threshold
        SWE_factor = SWEfactor_threshold; % if swe factor is abnormally high, adjust it to account for residue:
    end
    SWE_fbf_particles = SWE_pbp * SWE_factor; % adjusted SWE_pbp is effectivley SWE_fbf
    SWE_factor_particles = SWE_fbf_particles ./ SWE_pbp; % now calculate SWE factor for all particles
    SWE_pbp_accumulated = cumsum(SWE_pbp); % accumulated SWE_pbp
    SWE_fbf_accumulated = cumsum(SWE_fbf_particles); % accumulated SWE_fbf    
    snow_pbp = rho_water * (SWE_pbp ./ h_rho_hfd); % [mm]
    snow_pbp_accumulated = cumsum(snow_pbp); % [mm]
    snow_fbf = rho_water * (SWE_fbf_particles ./ h_rho_hfd); % [mm]
    snow_fbf_accumulated = cumsum(snow_fbf); % [mm]

    % store PBP data as table for post processing 
    pbp_table = table(h_initial_time', h_evap_time', ... 
        h_mass_pbp', h_diam',  h_max_area', h_ellipse_area', ... 
        h_rho_sph', h_rho_hfd', h_vol_hfd', h_vol_sph', ... 
        h_delta_temp_max', h_delta_temp_mean', complexity', ... 
        sdi', SWE_pbp', SWE_fbf_particles', SWE_pbp_accumulated', ... 
        SWE_fbf_accumulated', snow_pbp', snow_fbf', snow_pbp_accumulated', ... 
        snow_fbf_accumulated', SWE_factor_particles', deltaTemp_range', hArea_range', deltaTemp_residue_flags', hArea_residue_flags');
    pbp_table.Properties.VariableNames = {'Time', 'Evap Time (s)', 'Mass (kg)', 'Eff Diameter (m)',  ...
        'Max Area (m^2)', 'Ellipse Area (m^2)', 'Spherical Density (kg/m^3)', 'Heat Flux Density (kg/m^3)', ... 
        'Heat Flux Volume (m^3)', 'Spherical Volume (m^3)', 'Delta Temp Max', 'Delta Temp Mean', ...
        'Complexity','SDI','PBP SWE (mm)','FBF SWE (mm)', 'PBP SWE Accumulation (mm)', ...
        'FBF SWE Accumulation (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'PBP Snow Accumulation (mm)', ...
        'FBF Snow Accumulation (mm)', 'SWE factor', 'Delta Temp Range', 'Area Range', 'Delta Temp Flag', 'Area Flag'};
    pbp_table = table2timetable(pbp_table); % convert to timetable 
    pbp_table = sortrows(pbp_table, 'Time'); % sort by time 
    pbp_table.("PBP SWE Accumulation (mm)") = cumsum(pbp_table.("PBP SWE (mm)"));
    pbp_table.("FBF SWE Accumulation (mm)") = cumsum(pbp_table.("FBF SWE (mm)"));
    pbp_table.("PBP Snow Accumulation (mm)") = cumsum(pbp_table.("PBP Snow (mm)"));
    pbp_table.("FBF Snow Accumulation (mm)") = cumsum(pbp_table.("FBF Snow (mm)")); 
    pbp_table.("Missing Data") = false(height(pbp_table),1); % add a flag column to distinguish missing .avi data

    % handles FBF SWE data (used to compute SWE-factor later on)
    % fbf_table_raw = table(time_series_fbf', SWE_fbf);
    % fbf_table_raw.Properties.VariableNames{'Var1'} = 'Time';
    % fbf_table_raw.Time = datetime(fbf_table_raw.Time);
    % fbf_table_raw = table2timetable(fbf_table_raw);

    %% filtering starts here
    % filters data to find where 0 < mass < .005 and 1/15 s < evaporation time < 60 s to omit residue on plate
    f1 = find((pbp_table.("Mass (kg)") > 0 & ...
        pbp_table.("Mass (kg)") < residue_filter) & ...
        (pbp_table.("Evap Time (s)") > evapTime_min & ...
        pbp_table.("Evap Time (s)") < evapTime_max) & ...
        pbp_table.('Delta Temp Flag') ~= 1);

    pbp_table_filtered = pbp_table(f1,:);
    pbp_table_filtered.("PBP SWE Accumulation (mm)") = cumsum(pbp_table_filtered.("PBP SWE (mm)"));
    pbp_table_filtered.("FBF SWE Accumulation (mm)") = cumsum(pbp_table_filtered.("FBF SWE (mm)"));
    pbp_table_filtered.("PBP Snow Accumulation (mm)") = cumsum(pbp_table_filtered.("PBP Snow (mm)"));
    pbp_table_filtered.("FBF Snow Accumulation (mm)") = cumsum(pbp_table_filtered.("FBF Snow (mm)"));
    
    if height(pbp_table_filtered) > 0
        %% average last two minutes of data to account for time gap between .avi files 
        % grab last two minutes of data:
        final_time = pbp_table.Time(end);
        prev_time = final_time - minutes(2);
        prev_data = pbp_table(pbp_table.Time >= prev_time, :);
    
        % take the mean or sum of each value corresponding to the captured data:
        timeRow = final_time + seconds(5); 
        evapTimeRow = mean(prev_data.("Evap Time (s)"));
        massRow = sum(prev_data.('Mass (kg)'));
        effDiaRow = mean(prev_data.('Eff Diameter (m)'));
        areaRow =  mean(prev_data.('Max Area (m^2)'));
        ellipseAreaRow = mean(prev_data.("Ellipse Area (m^2)")); 
        volumeHFDrow = sum(prev_data.("Heat Flux Volume (m^3)"));
        volumeSPHrow = sum(prev_data.("Spherical Volume (m^3)"));
        deltaTempMaxRow = mean(prev_data.("Delta Temp Max"));
        deltaTempMeanRow = mean(prev_data.("Delta Temp Mean"));
        cxRow = mean(prev_data.Complexity);
        sdiRow = mean(prev_data.SDI);
        SWEfactorRow = mean(prev_data.("SWE factor"));
        tempRangeRow = NaN;
        aRangeRow = NaN; 
        tempFlagRow = NaN;
        aFlagRow = NaN; 
        densityHFDrow = massRow / volumeHFDrow;
        densitySPHrow = massRow / volumeSPHrow;
        PBPsweRow = 1000 * massRow/ (rho_water * hp_area);
        FBFsweRow = PBPsweRow*SWEfactorRow;
        PBPsnowRow = rho_water * (PBPsweRow ./ densityHFDrow);
        FBFsnowRow = rho_water * (FBFsweRow ./ densityHFDrow); 
        % compile into a timetable
        new_row = table(timeRow,... 
            evapTimeRow, massRow, effDiaRow, areaRow, ellipseAreaRow, ...
            densitySPHrow, densityHFDrow, volumeHFDrow, volumeSPHrow, deltaTempMaxRow, ...
            deltaTempMeanRow, cxRow, sdiRow, PBPsweRow, FBFsweRow, ...
            sum(pbp_table.("PBP SWE (mm)"))+PBPsweRow,...
            sum(pbp_table.("FBF SWE (mm)"))+FBFsweRow, ...
            PBPsnowRow, FBFsnowRow, ...
            sum(pbp_table.("PBP Snow (mm)"))+PBPsnowRow, ...
            sum(pbp_table.("FBF Snow (mm)"))+FBFsnowRow, SWEfactorRow, ...
            tempRangeRow, aRangeRow, tempFlagRow, aFlagRow, ...
            'VariableNames', {'Time', 'Evap Time (s)', 'Mass (kg)', 'Eff Diameter (m)',  ...
            'Max Area (m^2)', 'Ellipse Area (m^2)', 'Spherical Density (kg/m^3)', 'Heat Flux Density (kg/m^3)', ... 
            'Heat Flux Volume (m^3)', 'Spherical Volume (m^3)', 'Delta Temp Max', 'Delta Temp Mean', ...
            'Complexity','SDI','PBP SWE (mm)','FBF SWE (mm)', 'PBP SWE Accumulation (mm)', ...
            'FBF SWE Accumulation (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'PBP Snow Accumulation (mm)', ...
            'FBF Snow Accumulation (mm)', 'SWE factor',  'Delta Temp Range', 'Area Range', 'Delta Temp Flag', 'Area Flag'});
        % assign a time and logical value for missing data to the new row:
        new_row.('Missing Data') = true;
        % convert to timetable:
        new_row = table2timetable(new_row); 
        % now append new row to particle_output_table:
        pbp_table = [pbp_table; new_row];

        %%  create a summary table with actual data 
        avi_summary_table = table(pbp_table_filtered.Time(1));
        avi_summary_table.duration = (pbp_table_filtered.Time(end)-pbp_table_filtered.Time(1)); 
        avi_summary_table.cx = mean(pbp_table_filtered.Complexity);
        avi_summary_table.sdi = mean(pbp_table_filtered.SDI);
        avi_summary_table.rho = sum(pbp_table_filtered.("Mass (kg)")) / sum(pbp_table_filtered.("Heat Flux Volume (m^3)"));
        avi_summary_table.pbpSWE = pbp_table_filtered.("PBP SWE Accumulation (mm)")(end);
        avi_summary_table.fbfSWE = pbp_table_filtered.("FBF SWE Accumulation (mm)")(end);
        avi_summary_table.pbpSnow = pbp_table_filtered.("PBP Snow Accumulation (mm)")(end);
        avi_summary_table.fbfSnow = pbp_table_filtered.("FBF Snow Accumulation (mm)")(end);
        avi_summary_table.hotPlateArea = hp_area; 
        avi_summary_table.SWEfactor = SWE_factor; 
        avi_summary_table.minMass = h_mass_fbf_min;
    else 
        % create a summary table with all 0's
        avi_summary_table = table(pbp_table.Time(1));
        avi_summary_table.duration = (pbp_table.Time(end)-pbp_table.Time(1)); 
        avi_summary_table.cx = 0;
        avi_summary_table.sdi = 0;
        avi_summary_table.rho = 0;
        avi_summary_table.pbpSWE = 0;
        avi_summary_table.fbfSWE = 0;
        avi_summary_table.pbpSnow = 0;
        avi_summary_table.fbfSnow = 0;
        avi_summary_table.hotPlateArea = hp_area;
        avi_summary_table.SWEfactor = SWE_factor;
        avi_summary_table.minMass = h_mass_fbf_min;
    end
    % add variable names     
    % avi_summary_table = timetable2table(avi_summary_table);
    avi_summary_table.Properties.VariableNames = {'Time', 'Duration', 'Complexity', 'SDI', 'HFD Density (kg*m^-3)', 'PBP SWE (mm)', 'FBF SWE (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'Hot Plate Area', 'SWE Factor', 'Min FBF Mass'}; 
    avi_summary_table = table2timetable(avi_summary_table);

    % store as a cell so can be used outside of parforloop:
    pbp_table_cell{file_i} = pbp_table; 
    pbp_table_filtered_cell{file_i} = pbp_table_filtered;  
    avi_summary_table_cell{file_i} = avi_summary_table; 

 end

%% stored everything as a cell to access out of parfor loop; 
% now put back into one large time table 

% unfiltered particle data:
pbp_table = vertcat(pbp_table_cell{:});
pbp_table = sortrows(pbp_table, 'Time'); % sort by time
pbp_table.("PBP SWE Accumulation (mm)") = cumsum(pbp_table.("PBP SWE (mm)"));
pbp_table.("FBF SWE Accumulation (mm)") = cumsum(pbp_table.("FBF SWE (mm)"));
pbp_table.("PBP Snow Accumulation (mm)") = cumsum(pbp_table.("PBP Snow (mm)"));
pbp_table.("FBF Snow Accumulation (mm)") = cumsum(pbp_table.("FBF Snow (mm)")); 

% filtered particle data:
pbp_table_filtered = vertcat(pbp_table_filtered_cell{:});
pbp_table_filtered = sortrows(pbp_table_filtered, 'Time'); % sort by time
pbp_table_filtered.("PBP SWE Accumulation (mm)") = cumsum(pbp_table_filtered.("PBP SWE (mm)"));
pbp_table_filtered.("FBF SWE Accumulation (mm)") = cumsum(pbp_table_filtered.("FBF SWE (mm)"));
pbp_table_filtered.("PBP Snow Accumulation (mm)") = cumsum(pbp_table_filtered.("PBP Snow (mm)"));
pbp_table_filtered.("FBF Snow Accumulation (mm)") = cumsum(pbp_table_filtered.("FBF Snow (mm)")); 

% .avi summary table (uses filtered particle data):
avi_summary_table = vertcat(avi_summary_table_cell{:});
avi_summary_table = sortrows(avi_summary_table, 'Time'); % sort by time

% put .avi start times, end times, and lengths into one table:
% start_end_time_table = table(file_names', vid_start_time, vid_end_time); 
% start_end_time_table.length = start_end_time_table.vid_end_time - start_end_time_table.vid_start_time; 
% start_end_time_table.Properties.VariableNames{1} = 'file_name';

%% save processed tables

startTime = datestr(pbp_table.Time(1), 'yyyy-mm-dd_HH-MM-ss');

% unfiltered particle data table:
writetimetable(pbp_table, [output_dir, '/DEID_unfilteredParticle_v3_', startTime, '.csv']);

% filtered particle data table:
writetimetable(pbp_table_filtered, [output_dir, '/DEID_filteredParticle_v3_', startTime, '.csv']);

% .avi summary table:
writetimetable(avi_summary_table, [output_dir, '/DEID_aviTotals_v3_', storm_output, '.csv']);

% Writes out time averaged data table 
% writetimetable(ts_output_table, [output_dir,'/DEID_TS_10min_', startTime, '.csv']);

[~, parent_dir, ~] = fileparts(pwd);
disp(['Saved Output for: ', parent_dir])