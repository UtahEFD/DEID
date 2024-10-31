%% DEID AVI File Processing Code
% Outputs particle by particle csv file to be post processed by
% DEID_AVI_Processor.m
% AUTHOR : Dhiraj Singh, Benjamin Silberman, Travis Morrison, Alex Blackmer
% LAST UPDATED: 05/08/2024
                                                                 
clear, clc, close all
%% Sets filepath, global variables, and physical constants
working_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/JAN/JAN1';
output_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/Parsivel_DEID_Comparison/DEID_Data/v2/2min/';
% working_dir = 'Z:\DEID\Atwater\JAN\test';     % For use on Snowpack

% specifies resampling period:
time_interval = 120;  % Seconds
time_step = seconds(time_interval); % Datetime step 
% unit conversions:
mm_to_inches = 1/25.4; % [mm/in]
pix_to_m_conversion = 3.1750e-04; % Pixel to [m]
pix_to_m2_conversion = pix_to_m_conversion^2; % Pixel to [m^2]
k = 3.6e06; % Converts m/s to mm/hr [(mm/hr)*(s/m)]
hour_in_seconds = 3600;
% physical constants:
rho_water = 1000; % Density of water [kg/m^3]
mu = 1.5*10^-5;   % Viscosity of air [kg/m*s] 
% DEID specific parameters:
residue_filter = 0.005; % [kg]
colorbar_image_indexes = [1 1 384 288]; % Location of colorbar in pixel locations
colorbar_kapton_image_indexes = [1 39 383 261]; % Location of Kapton tape in pixel locations
colorbar_max_temp = 145; % Max temperature set in colorbar on the physical screen of the tir software
min_thres = 70; % Minimum threshold number in image accepted rbg ([0 255]) scale
sort_threshold = 20; % This it the RMS threshold between succesive images of snowflakes used in the sortPostitions_v2
min_h_size = 10; % Minumum Hydrometeor size in pixels 
minimum_drop_life = 0; % Minimum number of frames a drop has to be visable to be processed
k_dLv = 0.0035; % Calibration constant, in paper thermal conductivity (k) of water See sect 4.1 in Dihiraj's paper -> (k/d(_eff))/Latent heat of vaporazation [units?]
l_constant = 2.594e06; % Latent heat of vaporazation of water, should be a function of tempertaure (Look at Stull textbook) [J/kg]
% Eqn. (13) in Dhiraj's density paper: c = (L_vv) / (L_ff*C_melt)
% c = hf_rho_coeff
hf_rho_coeff = 1.01e05; % [K*s*m^-1] 

%% Move to working directory and identify video files 
% Sets path for .m to run in:
cd(working_dir) 
% Get a list of all files and folders in this directory
directory = dir(".");
% Initialize an empty cell array to store the file names
file_names = {};
% Loop through each item in the input directory
for file_i = 1:length(directory)
    % Check if the item is a file (not a folder) and if it ends with .avi
    if ~directory(file_i).isdir && endsWith(directory(file_i).name, '.avi', 'IgnoreCase', true)
        % Get the name of the file and append it to the list
        file_names{end+1} = directory(file_i).name;
    end
end

%% Now loop through each file_name to get start and end times 

vid_end_time = NaT(length(file_names), 1);
vid_start_time = NaT(length(file_names), 1);
vid_end_start_diff = NaN(length(file_names)-1, 1);

parfor file_i = 1:length(file_names)
    filename = file_names{file_i};
    vid=VideoReader(filename);
    % Get necessary metadata for video processing
    vid_dir = dir(filename);
    vid_length = vid.Duration;
    vid_fps = vid.FrameRate;
    num_frames = vid.NumFrames;
    % Create a time series of date times starting from (video end time - video duration) and ending at the video end time
    vid_end_time(file_i) = datetime([vid_dir.date]);
    time_series = datetime(vid_end_time(file_i) - (0:num_frames) * seconds(1/vid_fps), 'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
    time_series = flip(time_series);  % Flips time series to be chronologically ordered
    vid_start_time(file_i) = datetime(time_series(1));
end

start_end_time_table = table(file_names', vid_start_time, vid_end_time); 
start_end_time_table.length = start_end_time_table.vid_end_time - start_end_time_table.vid_start_time; 
start_end_time_table.Properties.VariableNames{1} = 'file_name'; 

% Find difference between end of one video and the start of the next
for file_i = 1:length(file_names)-1
    vid_end_start_diff(file_i) = seconds(start_end_time_table.vid_start_time(file_i+1) - start_end_time_table.vid_end_time(file_i));
end

%% Initialize Output Tables
% Frame by Frame output table 
% fbf_col_names = {'Time', 'SWE_mm'};
% fbf_col_types = {'datetime', 'double'};
% fbf_output_table = table('Size', [0, length(fbf_col_names)], ...
%                          'VariableNames', fbf_col_names, ...
%                          'VariableTypes', fbf_col_types);
% Time series output table
ts_col_names = {'Time', 'Terminal_Velocity', 'Complexity', 'SDI', 'Mass', 'Volume', 'Density_HFD', 'Diameter', 'Surface Area', 'Void Space', 'Temp Diff', 'SWE_mm','Snow_mm'};
ts_col_types = {'datetime', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double' ,'double', 'double'};
ts_output_table = table('Size', [0, length(ts_col_names)], ...
                         'VariableNames', ts_col_names, ...
                         'VariableTypes', ts_col_types);
% Particles output table
particle_col_names = {'Time', 'Terminal_Velocity', 'Complexity', 'SDI', 'Mass', 'Volume', 'Density_HFD', 'Diameter', 'Surface Area', 'Void Space', 'Temp Diff', 'SWE_mm', 'Snow_mm', 'Missing Data'};
particle_col_types = {'datetime', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double' ,'double', 'double', 'logical'};
particle_output_table = table('Size', [0, length(particle_col_names)], ...
                         'VariableNames', particle_col_names, ...
                         'VariableTypes', particle_col_types);
% Diagnostic output table
diag_col_names = {'Filename', 'Start_Time','End_Time', 'SWE_Factor', 'Num_Particles'};
diag_col_types = {'string', 'datetime', 'datetime', 'double', 'double'};
diag_output_table = table('Size', [0, length(diag_col_names)], ...
                         'VariableNames', diag_col_names, ...
                         'VariableTypes', diag_col_types);

%% Begin DEID video processing:

% % if processing a storm, specify storm start and end date:
% storm_start = datetime('07-Jan-2024 04:33:40');  
% storm_end = datetime('08-Jan-2024 07:26:38'); 
% storm_table = start_end_time_table(start_end_time_table.vid_start_time >= storm_start & start_end_time_table.vid_end_time <= storm_end, :); 
% file_names = storm_table.file_name; 
%%
% Parfor loop parallelizes processing by distributing each video file to a
% Matlab worker on each CPU core. 

final_particle_output_table = table(); 
hp_area = NaN(1, length(file_names)); % Preallocate Hotplate Area [m^2]
swe_factor = NaN(1, length(file_names)); % Preallocate SWE factor 

parfor file_i = 1:length(file_names)
    diag_table = table('Size', [1, length(diag_col_names)], ...
                         'VariableNames', diag_col_names, ...
                         'VariableTypes', diag_col_types);
    filename = file_names{file_i};
    disp(['Processing File: ', filename])
    vid=VideoReader(filename);
    
    % Creates datetime timeseries for output file
    % Get necessary metadata for video processing
    vid_dir = dir(filename);
    vid_length = vid.Duration;
    vid_fps = vid.FrameRate;
    num_frames = vid.NumFrames;

    % Gets video start and end date to construct time series
    vid_end_time = datetime([vid_dir.date]);
    
    % Create a time series of date times starting from (video end time - video duration) and ending at the video end time
    time_series = datetime(vid_end_time - (0:num_frames) * seconds(1/vid_fps), 'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
    time_series = flip(time_series);  % Flips time series to be chronologically ordered
    vid_start_time = datetime(time_series(1));

    % start_end_time_table = tablength(file_names)le(vid_start_time, vid_end_time)
    % Get number of frames, preallocate cell for data, and obtain date information:
    h_data = cell(num_frames,1);    % Hydrometeor dat

    % Appends diagnostic data to table
    diag_table.Filename = filename;
    diag_table.Start_Time = vid_start_time;
    diag_table.End_Time = vid_end_time;

    %% "Frame by frame method"; this is how Dhiraj is obtaining SWE for each .avi file 
    % Preallocate variables saved in loop for speed:
    plate_temp = nan(num_frames,1);
    sum_h_area_times_dt = nan(num_frames,1);
    % Enter loop to process images: 
    for frame_ii = 1:num_frames 
        % Get image:    
        frame = read(vid, frame_ii); 
        % Clean image to get plate temperature: 
        frame_gray = im2gray(frame); % Convert frame of interest to gray scale
        frame_gray_cropped_wKapton = imcrop(frame_gray, colorbar_image_indexes);% Crop out colorbar
        plate_temp(frame_ii) = max(max(double(frame_gray_cropped_wKapton))); % Dhiraj assumes max temperature in image is the plate temperature with Kapton tape
        % Clean orginal image: 
        frame_gray_cropped = imcrop(frame_gray, colorbar_kapton_image_indexes); % Back to orginal grayscale image... now remove colorbar and kapton tape from image
        frame_filtered = frame_gray_cropped > min_thres; % Removed below min threshold, on rbg ([0, 255]) scale 
        frame_filtered_filled = imfill(frame_filtered, 'Holes'); % Clean up Hydrometeors
        % Get Hydrometeor properties: 
        h_geo_prop = regionprops(frame_filtered_filled, 'Centroid', 'Area','BoundingBox'); % Returns the centroid, the area of each blob, and the bounding box (left, top, width, height).
        % If no properties are found, go to next frame: 
        if (isempty(h_geo_prop))
            continue;
        end
        % Build Hydrometeor propertie matrices from regionprops values: 
        h_bounding_box = cat(1,h_geo_prop.BoundingBox); % Concat all values to build matrix of Bounding Box indices
        h_centroid = round(cat(1, h_geo_prop.Centroid)); % Concat all values to build matrix of Centroids
        h_area_pix = cat(1, h_geo_prop.Area); % Concat all values to build matrix of Hydrometeor areas in pixels
        % Convert Hydrometeor areas to m^2: 
        h_area = h_area_pix .* pix_to_m2_conversion; 
        h_major_axis = h_bounding_box(:, 3) * sqrt(pix_to_m2_conversion); % Hydrometeor major axis
        h_minor_axis = h_bounding_box(:, 4) * sqrt(pix_to_m2_conversion); % Hydrometeor minor axis
        % Get the Hydrometeor ellipse area:
        h_elipse_area = h_major_axis .* h_minor_axis;
        % Calculating the difference in temperature of each centroid and
        % the plate:
        h_centroid_i = sub2ind(size(frame_gray_cropped), h_centroid(:, 2), h_centroid(:, 1)); % Find the index of the centriods in orginal image
        h_centroid_values = double(frame_gray_cropped(h_centroid_i)); % Intensities of Centroid pixels of snow
        plate_h_dtemp = colorbar_max_temp - (h_centroid_values .* (colorbar_max_temp / plate_temp(frame_ii)));
        % Now calculate product of Hydrometeor area with the temp
        % difference:
        h_area_times_dtemp = h_area .* plate_h_dtemp;         
        % Sum the product of individual area and temp. diff in each frame:
        % THIS IS HOW WE GET Hydrometeor MASS (FRAME BY FRAME METHOD)
        sum_h_area_times_dt(frame_ii)=sum(h_area_times_dtemp);         
        % Build large matrix of Hydrometeor data
        h_data{frame_ii} = cat(2, h_centroid, h_area, plate_h_dtemp, h_elipse_area, h_major_axis, h_minor_axis); 
    end

    % Frame by frame SWE calculation
    % Use current_frame_gray_cropped to calculate the area of hot plate:
    hp_area(file_i) = size(frame_gray_cropped,1) * size(frame_gray_cropped,2) * pix_to_m2_conversion;    
    % hp_area(file_i) = hp_area_temp; % Hotplate Area [m^2]
    h_mass_fbf = (k_dLv*sum_h_area_times_dt) / vid_fps; % mass evaporates in each frame
    SWE_FBF_mm = h_mass_fbf / hp_area(file_i);
    % Find the minimum SWE in all frames within a video, and subtract from
    % SWE (way of handling residue) 
    % SWE_FBF_mm = SWE_FBF_mm - min(SWE_FBF_mm);
    SWE_fbf_accumulation = sum(SWE_FBF_mm); 
    time_series_fbf = time_series(1:length(SWE_FBF_mm));
    
    %% Sorting one frame to others frame ~ Data cleaning of some sorts: 
    % Create new cell array for sorted data: 
    h_data_sorted = cell(size(h_data));
    h_data_sorted{1} = h_data{1};
    % Loop over every data cell corresponding to each frame:
    % Uses "sortPositions_v2.m" - a code that Dhiraj wrote 
    for frame_jj = 2:num_frames
        h_data_sorted{frame_jj} = h_data{frame_jj};
        h_data_sorted{frame_jj} = sortPositions_v2(h_data_sorted{frame_jj-1}, h_data_sorted{frame_jj}, sort_threshold);
    end

    % Return frame with max number of Hydrometeors 
    max_h_obs = max(cellfun(@(x) size(x, 1), h_data_sorted, 'UniformOutput', 1));
    % Now pad the data with zeros so the frames all have the same number:  
    h_data_sorted = cellfun(@(x) cat(1, x, zeros(max_h_obs - size(x, 1), 7)),...
        h_data_sorted, 'UniformOutput', 0);
    
    %% Isolating the variables and put them into a matrix to work with:
    % For reference: [Hydrometeor_Centroid, Hydrometeor_Area,
    % Plate_Hydrometeor_DeltaT, Hydrometeor_ellipse_area,
    % Hydrometeor_Major_Axis,Hydrometeor_Minor_Axis]
    h_area_final = cellfun(@(x) x(:, 3), h_data_sorted, 'UniformOutput', 0);
    dT_final = cellfun(@(x) x(:, 4), h_data_sorted, 'UniformOutput', 0);
    ellipse_area_final = cellfun(@(x) x(:, 5), h_data_sorted, 'UniformOutput', 0);
    maj_axis_final = cellfun(@(x) x(:, 6), h_data_sorted, 'UniformOutput', 0);
    min_axis_final = cellfun(@(x) x(:, 7), h_data_sorted, 'UniformOutput', 0);
    
    % Convert to matrix: 
    h_area_final = cat(2,h_area_final{:});
    dT_final = cat(2,dT_final{:});
    ellipse_area_final = cat(2,ellipse_area_final{:});
    maj_axis_final = cat(2,maj_axis_final{:});
    min_axis_final = cat(2,min_axis_final{:});

    %% "Particle by particle method"
    % Calulcate the evolution of each Hydrometeor to get properties.
    % This method 'cleans' the plate to not double count Hydrometeors that
    % appear across multiple frames. 
    % Preallocate variables:
    h_delta_time = {}; % Hydrometeor time to melt and evaporate
    h_delta_temp_max = {}; % Hydrometeor centroid's max intensity
    h_delta_temp_mean = {}; % Hydrometeor centroid's mean intensity
    h_max_area = {}; % Hydrometeor maximum area
    h_max_circ_area = {};% Hydrometeor Max circumscribed area
    h_max_maj_axis = {};% Hydrometeor Max major axis
    h_max_min_axis = {};% Hydrometeor Max minor axis
    h_mass_pbp = {}; % Hydrometeor calculated mass
    h_density = {}; % Hydrometeor density
    h_init_time_ind = []; % Initial time index of each snowflake in time array

    % Loop over all Hydrometeors: 
    % This is where the cleaning occurs!
    for h_ii = 1:max_h_obs
        h_appear_evap_bool = diff(h_area_final(h_ii, :) > 0); 
        % Create array of indexes which are booleans, +1 is when snowflake starts, and -1 when it leaves
        % diff is [X(3)-X(2)] for backwward FD
        
        h_appears_ind = find(h_appear_evap_bool > 0); % finds the index when the Hydrometeor appears
        if isempty(h_appears_ind)
            continue;
        end
        h_evaps_ind = find(h_appear_evap_bool < 0); % Hydrometeor dissapears index
        h_init_time_ind = cat(2,h_init_time_ind,h_appears_ind); % time when snow arrives the plate
        if isempty(h_evaps_ind) % get rid of Hydrometeors which do no evaporate completely
            continue; 
        else
        % Isolate the Hydrometeor properties for the lifetime of the
        % Hydrometeor: 
        h_dT_tmp = dT_final(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
        h_ellipse_area_tmp = ellipse_area_final(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
        h_maj_axis_tmp = maj_axis_final(h_ii,h_appears_ind(1):h_evaps_ind(end)+1);
        h_min_axis_tmp = min_axis_final(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
        % Prepare area of each Hydrometeor for temporal intergration: 
        h_area_tmp = h_area_final(h_ii, h_appears_ind(1):h_evaps_ind(end)+1); 
        h_area_tmp_bool = h_area_tmp > 0; % Isolate for positive values
        h_area_tmp_bool = bwareaopen(h_area_tmp_bool, minimum_drop_life); % Any Hydrometeor which lives for less than 'Minimum_Drop_Life' is discarded
        % Integration ~ need to double check with google search: 
        propstemp = regionprops(h_area_tmp_bool, 'PixelIdxList');
        % Loop over some region property definded in previous line: 
            for jj = 1:numel(propstemp)
                h_max_area{end+1} = max(h_area_tmp(propstemp(jj).PixelIdxList));
                h_max_circ_area{end+1} = max(h_ellipse_area_tmp(propstemp(jj).PixelIdxList));
                h_max_maj_axis{end+1} = max(h_maj_axis_tmp(propstemp(jj).PixelIdxList));
                h_max_min_axis{end+1} = max(h_min_axis_tmp(propstemp(jj).PixelIdxList));
                h_delta_time{end+1} = numel(propstemp(jj).PixelIdxList);
                h_delta_temp_max{end+1} = max(h_dT_tmp(propstemp(jj).PixelIdxList));
                h_delta_temp_mean{end+1} = mean(h_dT_tmp(propstemp(jj).PixelIdxList));
                % Multiply the area by the temperature difference for that Hydrometeor:
                h_area_times_dT_pbp = h_area_tmp(propstemp(jj).PixelIdxList) .* h_dT_tmp(propstemp(jj).PixelIdxList);
                % Now integrate:
                h_mass_pbp{end+1} = sum(k_dLv * h_area_times_dT_pbp);
            end
        end
    end
    
    %% Convert to matrixes:
    h_mass_pbp=cell2mat(h_mass_pbp); % Hydrometeor mass
    h_max_area=cell2mat(h_max_area); % Actual area
    h_max_circ_area=cell2mat(h_max_circ_area); % Ellipse area
    h_max_maj_axis=cell2mat(h_max_maj_axis); % Height
    h_max_min_axis=cell2mat(h_max_min_axis); % Width
    h_delta_temp_max=cell2mat(h_delta_temp_max); % Temperature diffrence between plate and water droplet using max intensity
    h_delta_temp_mean=cell2mat(h_delta_temp_mean); % Temperature diffrence between plate and water droplet using mean intensity 
    % Now multiples mass by "dt":
    h_mass_pbp = h_mass_pbp / vid_fps; 
    h_diam = (4/pi) * h_max_area.^(1/2); % Convert area to diameter of Hydrometeor
    h_water_eq_diameter = (6 .* h_mass_pbp / (pi * rho_water)).^(1/3); % Water equi. diameter. Not sure where this is from
    h_evap_time = cell2mat(h_delta_time) * (1 / vid_fps); % Evaporation time
    h_sph_vol = (3/4) * h_max_area.^(3/2); % Spherical volume
    h_rho_sph = h_mass_pbp ./ h_sph_vol; % Density calculation: spherical assumption
    h_energy_per_time = h_mass_pbp * l_constant ./ (h_max_area .* h_evap_time); % Heat flux method: energy per unit area per time
    h_rho_hfd = (hf_rho_coeff * h_mass_pbp) ./ (h_max_area .* h_evap_time .* h_delta_temp_mean); 
    h_vol_hfd = h_mass_pbp ./ h_rho_hfd; % Volume of each snowflakes using mean heat flux method density
    h_initial_time_indexes = h_init_time_ind(1:length(h_mass_pbp));
    h_initial_time = time_series(h_initial_time_indexes);

    %% Organizes data into tables used for both output and post processing
    % Handles FBF SWE data
    fbf_table_raw = table(time_series_fbf', SWE_FBF_mm);
    fbf_table_raw.Properties.VariableNames{'Var1'} = 'Time';
    fbf_table_raw.Time = datetime(fbf_table_raw.Time);
    fbf_table_raw = table2timetable(fbf_table_raw);

    % Handles PBP data
    pbp_table_particles = table(h_initial_time', h_evap_time', h_mass_pbp', h_diam',  ...
        h_max_area', h_max_circ_area', h_rho_sph', h_rho_hfd',  ... 
        h_vol_hfd', h_water_eq_diameter', h_max_maj_axis', h_max_min_axis', ...
        h_energy_per_time', h_delta_temp_max', h_delta_temp_mean');
    pbp_table_particles.Properties.VariableNames = {'initial_time', 'evap_time', 'mass', 'diam',  ...
        'max_area', 'max_circ_area', 'rho_sph', 'rho_hfd', ... 
        'vol_hfd', 'water_eq_diameter', 'max_maj_axis', 'max_min_axis', ...
        'energy_per_time', 'delta_temp_max', 'delta_temp_mean'};
    
    pbp_table_particles = table2timetable(pbp_table_particles); % Convert to timetable 
    pbp_table_particles = sortrows(pbp_table_particles, 'initial_time'); % Sort by time 
    diag_table.Num_Particles = size(pbp_table_particles,1);

    %% Post Processing Starts Here! 
    % Processes PBP data if particles were present in video file
    if ~isempty(pbp_table_particles)

        % Filters data to find where 0 < mass < .005 to omit residue on plate
        % Change this to filter on evap time
        [g1,g2] = find(pbp_table_particles.mass > 0 & pbp_table_particles.mass < residue_filter); 
        pbp_table_particles = pbp_table_particles(g1,:);

        % Calculate terminal velocity using terminalVelocity function
        pbp_table_particles.terminal_vel = terminalVelocity(pbp_table_particles.max_area, pbp_table_particles.max_circ_area, pbp_table_particles.mass);
        
        % Calculate 'void space'
        pbp_table_particles.void_space = (pbp_table_particles.max_circ_area - pbp_table_particles.max_area).^(1/2); 

        % Calculate SDI, Cx  
        pbp_table_particles.surface_area_eq = ((9*pi)/16)^(1/3) * (pbp_table_particles.mass ./ rho_water).^(2/3); % Surface area of equivalent water droplet (See POF Singh et al. 2023)
        pbp_table_particles.complexity = pbp_table_particles.max_circ_area ./ pbp_table_particles.max_area; % 'Complexity' (See CRST Morrison et al. 2023)
        pbp_table_particles.sdi = pbp_table_particles.max_area ./ pbp_table_particles.surface_area_eq ; % 'SDI' (See CRST Morrison et al. 2023)

        % Resamples time series at desired interval
        fbf_table_retimed = retime(fbf_table_raw, 'regular', 'sum', 'TimeStep', time_step);
        fbf_table_retimed.SWE_FBF_accum_mm = cumsum(fbf_table_retimed.SWE_FBF_mm);

        %% Total SWE for all PBP data
        pbp_table_particles.SWE_PBP_mm = 1000 * pbp_table_particles.mass ./ (rho_water * hp_area(file_i)); % [mm]
        pbp_table_particles.SWE_PBP_accum_mm = cumsum(pbp_table_particles.SWE_PBP_mm);  % [mm]
     
        % Finds difference factor between FBF SWE and PBP SWE and adjusts PBP SWE
        swe_factor(file_i) = fbf_table_retimed.SWE_FBF_accum_mm(end) / pbp_table_particles.SWE_PBP_accum_mm(end);

        % Adjusts PBP SWE
        pbp_table_particles.SWE_PBP_F_mm = pbp_table_particles.SWE_PBP_mm * swe_factor(file_i); 
        pbp_table_particles.SWE_PBP_F_accum_mm = cumsum(pbp_table_particles.SWE_PBP_F_mm);
        
        %% Total Snow for all PBP data
        pbp_table_particles.snow_PBP_mm = rho_water * (pbp_table_particles.SWE_PBP_F_mm ./ pbp_table_particles.rho_hfd); % [mm]
        pbp_table_particles.snow_PBP_acc_mm = cumsum(pbp_table_particles.snow_PBP_mm); % [mm]

        %% Appends PARTICLE data for single video to output table
        % selects a subset of output variables to be exported:
        particle_output_table_all = pbp_table_particles(:, {'terminal_vel', 'complexity', 'sdi', 'mass', 'vol_hfd', 'rho_hfd', 'diam', 'surface_area_eq', 'void_space', 'delta_temp_mean', 'SWE_PBP_F_mm', 'snow_PBP_mm'});
        particle_output_table_all = timetable2table(particle_output_table_all);
        % add a flag column to distinguish missing .avi data:
        particle_output_table_all.missing_data = false(height(particle_output_table_all),1); 
        particle_output_table_all.Properties.VariableNames = particle_col_names;
        
        % Loop through each variable and replace NaNs with zeros.
        % NaNs are showing up for some videos because the pbp and fbf tables
        % dont have the same number of time stamps. Replacing Nans is a quick
        % fix but should be further explored.
        % particle by particle
        output_names = particle_output_table_all.Properties.VariableNames;
        for i = 2:length(output_names)
            particle_output_table_all.(output_names{i})(isnan(particle_output_table_all.(output_names{i}))) = 0;
        end

         % Appends output
        particle_output_table = [particle_output_table; particle_output_table_all(:, particle_col_names)];
        
        %% Find an average of previous two minutes of data to fill time gap between .avi files 

        % Go back two minutes and capture corresponding data:
        final_time = particle_output_table_all.Time(end);
        prev_time = final_time - minutes(2);
        prev_data = particle_output_table_all(particle_output_table_all.Time >= prev_time, :);

        % take the mean or sum of each value corresponding to the captured data:
        termVelocityRow = mean(prev_data.Terminal_Velocity);
        cxRow = mean(prev_data.Complexity);
        sdiRow = mean(prev_data.SDI);
        massRow = mean(prev_data.Mass);
        volumeRow = mean(prev_data.Volume);
        densityRow = massRow / volumeRow; 
        diamRow = mean(prev_data.Diameter);
        surAreaRow = mean(prev_data.('Surface Area'));
        voidSpaceRow = mean(prev_data.('Void Space'));
        tempDiffRow = mean(prev_data.('Temp Diff'));
        sweRow = sum(prev_data.('SWE_mm'));
        snowRow = sum(prev_data.('Snow_mm'));
        % compile into a table
        new_row = array2table([termVelocityRow, cxRow, sdiRow, massRow, ...
            volumeRow, densityRow, diamRow, surAreaRow, voidSpaceRow, tempDiffRow, ...
            sweRow, snowRow], 'VariableNames', particle_output_table_all.Properties.VariableNames(2:(end-1)));
        % assign a time and logical value for missing data to the new row:
        new_row.('Missing Data') = true;
        new_row.Time = final_time + seconds(5);
        % if file_i  == 1
        %     new_row = new_row(:, [end, 1:end-1]);
        % end

        % now append new row to particle_output_table:
        particle_output_table = [particle_output_table; new_row];

    end
end

%% Call missingTimeFunction to fill gaps between .avi files

% [full_storm_data_table] = missingTimeFunction(particle_output_table);
%% Cumulatively sums data for SWE and Snow totals
particle_output_table.SWE_Accum_mm = cumsum(particle_output_table.SWE_mm);
particle_output_table.Snow_Accum_mm = cumsum(particle_output_table.Snow_mm);

% full_storm_data_table.SWE_Accum_mm = cumsum(full_storm_data_table.SWE_mm);
% full_storm_data_table.Snow_Accum_mm = cumsum(full_storm_data_table.Snow_mm);

%% Averaging technique(s) starts here! 

particle_output_table = table2timetable(particle_output_table);
%% Computes averaged and summed PBP data using retime function
        
% Averages PBP data 
avg_cols = {'Terminal_Velocity', 'Complexity', 'SDI', 'Diameter', 'Surface Area', 'Void Space', 'Temp Diff'};
avg_table = retime(particle_output_table(:, avg_cols), 'regular', 'mean', 'TimeStep', time_step);
% Sums PBP data
sum_cols = {'Mass', 'Volume'};
sum_table = retime(particle_output_table(:, sum_cols), 'regular', 'sum', 'TimeStep', time_step);
% Join tables
pbp_table_retimed = horzcat(avg_table, sum_table);
        
%% Total SWE per averaging period PBP data
pbp_table_retimed.SWE_PBP_mm = 1000 * pbp_table_retimed.Mass ./ (rho_water * hp_area(1)); % [mm]
pbp_table_retimed.SWE_PBP_accum_mm = cumsum(pbp_table_retimed.SWE_PBP_mm);  % [mm]

% Finds difference factor between FBF SWE and PBP SWE and adjusts PBP SWE
% swe_factor = fbf_table_retimed.SWE_FBF_accum_mm(end) / pbp_table_retimed.SWE_PBP_accum_mm(end);
        
% Adjusts averaged PBP SWE
pbp_table_retimed.SWE_PBP_F_mm = pbp_table_retimed.SWE_PBP_mm * swe_factor(1);                                   
pbp_table_retimed.SWE_PBP_F_accum_mm = cumsum(pbp_table_retimed.SWE_PBP_F_mm);                   
        
% diag_table.SWE_Factor = swe_factor     

%% Total Snow per averaging period PBP data
% Adjusted density
pbp_table_retimed.rho_hfd = pbp_table_retimed.Mass ./ pbp_table_retimed.Volume;   % Density from HFD density method [kg/m^3]
pbp_table_retimed.snow_PBP_mm = rho_water * (pbp_table_retimed.SWE_PBP_mm ./ pbp_table_retimed.rho_hfd); % [mm]
pbp_table_retimed.snow_PBP_acc_mm = cumsum(pbp_table_retimed.snow_PBP_mm); % [mm]
        
%% Appends AVERAGED TIME SERIES data for single video to output table
% video_output_table = synchronize(fbf_table_retimed, pbp_table_retimed);
% Selects a subset of output variables to be exported
video_output_table = pbp_table_retimed(:, {'Terminal_Velocity', 'Complexity', 'SDI', 'Mass', 'Volume', 'rho_hfd', 'Diameter', 'Surface Area', 'Void Space', 'Temp Diff', 'SWE_PBP_F_mm','snow_PBP_mm'});
video_output_table = timetable2table(video_output_table);
video_output_table.Properties.VariableNames = ts_col_names;
        
% Loop through each variable and replace NaNs with zeros.
% NaNs are showing up for some videos because the pbp and fbf tables
% dont have the same number of time stamps. Replacing Nans is a quick
% fix but should be further explored.
output_names = video_output_table.Properties.VariableNames;
for i = 2:length(output_names)
    video_output_table.(output_names{i})(isnan(video_output_table.(output_names{i}))) = 0;
end
        
% Appends output
ts_output_table = [ts_output_table; video_output_table(:, ts_col_names)];
% diag_output_table = [diag_output_table; diag_table(:, diag_col_names)];

%% Sorts table by time and handles duplicates 
% particle_output_table = table2timetable(sortrows(particle_output_table, 'Time'));
ts_output_table = sortrows(ts_output_table, 'Time');

% Custom function averages some variables and sums others
customFunction = @(x) [mean(x(:,1:10), 1), sum(x(:,11:12), 1)]; 
averagedAndSummedValues = splitapply(customFunction, table2array(ts_output_table(:,2:end)), ...
                                     findgroups(ts_output_table.Time));
ts_output_table = table2timetable(array2table(averagedAndSummedValues, ...
                                            'VariableNames', ts_output_table.Properties.VariableNames(2:end)), ...
                               'RowTimes', unique(ts_output_table.Time));

% diag_output_table = sortrows(diag_output_table, 'Start_Time');
%% Cumulatively sums data for SWE and Snow totals
% Averaged data 
ts_output_table.SWE_Accum_mm = cumsum(ts_output_table.SWE_mm);
ts_output_table.Snow_Accum_mm = cumsum(ts_output_table.Snow_mm);
ts_output_table.Snow_Accum_in = ts_output_table.Snow_Accum_mm * mm_to_inches;

%% Now create an averaged data table using a given snow interval

%% Saves processed output and diagnostic data for all video files present
% Gets folder name and saves output as 'folder name'.csv
startTime = datestr(ts_output_table.Time(1), 'yyyy-mm-dd_HH-MM-ss');

% Writes out all particle data table
writetimetable(particle_output_table, [output_dir, 'DEID_Particle_2min_', startTime, '.csv']);

% Writes out particle data table including missing time between .avi files
% writetimetable(full_storm_data_table, ['DEID_missingParticle_TEST_', startTime, '.csv']);

% Writes out time averaged data table 
writetimetable(ts_output_table, [output_dir,'DEID_TS_2min_', startTime, '.csv']);

% Writes out snow interval ageraged data table
% writetimetable(snowInterval_table, ['DEID_snowAvg_', startTime, '.csv']);
 
% Writes out diagnostic data table
writetable(diag_output_table, [output_dir,'Diag_DEID_', startTime,'.csv']);

[~, parent_dir, ~] = fileparts(pwd);
disp(['Saved Output for: ', parent_dir])