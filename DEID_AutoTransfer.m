%% DEID AVI File Processing Code
% Reads in individual .avi files to push to chpc & website 
% AUTHOR : Dhiraj Singh, Benjamin Silberman, Travis Morrison, Alex Blackmer

clear, clc
%% Sets filepath, global variables, and physical constants.
working_dir = 'D:\Atwater';
% working_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/JAN/JAN1';     % For testing
output_dir = 'D:\Atwater\feb19';

% Set global varables and constants:
% specifies resampling period:
time_interval = 300;  % Seconds
time_step = seconds(time_interval); % Datetime step 
% unit conversions:
mm_to_inches = 1/25.4; % [mm/in]

% pix to m conversion:
hotPlateArea_m2 = 0.0130; % measuring from IRBIS 3.1 plus software 
hotPlateArea_pixels = 110592.00; % measuring from IRBIS 3.1 plus software
pix_to_m2_conversion = hotPlateArea_m2/hotPlateArea_pixels;
c1 = 10^3; % Converts meters to millimeters
% pix_to_m_conversion = 3.1750e-04; % pixel to [m]
% pix_to_m2_conversion_old = pix_to_m_conversion^2; % pixel to [m^2]
k = 3.6e06; % converts m/s to mm/hr [(mm/hr)*(s/m)]
hour_in_seconds = 3600;
% physical constants:
rho_water = 1000; % density of water [kg/m^3]
mu = 1.5*10^-5;   % viscosity of air [kg/m*s] 
l_v = 2.32e06; % latent heat of vaporazation of water [J/kg]
l_vv = 2.594e06; % what is this? l_vv? 
l_f = 3.34e05; % latent heat of fusion of water [J/kg]
% DEID specific parameters:
residue_filter = 0.005; % [kg]
evapTime_min = 1/15;
% evapTime_save = '_one'; 
evapTime_max = 60;
colorbar_image_indexes = [1 1 384 288]; % Location of colorbar in pixel locations
crop_index = 43; % use this to specify indices to crop out kapton tape
colorbar_kapton_image_indexes = [1 (colorbar_image_indexes(2)+crop_index) 383 (colorbar_image_indexes(4)-crop_index)]; % Location of Kapton tape in pixel locations
% [x1 y1 width height] y1 + n, height - n
colorbar_max_temp = 145; % Max temperature set in colorbar on the physical screen of the tir software
min_thres = 70; % Minimum threshold number in image accepted rbg ([0 255]) scale
sort_threshold = 20; % This it the RMS threshold between succesive images of snowflakes used in the sortPostitions_v2
min_h_size = 10; % Minumum Hydrometeor size in pixels 
minimum_drop_life = 0; % Minimum number of frames a drop has to be visable to be processed
k_dEff = 7.006e3; % calibration constant (see sect 4.1 of dihiraj's DEID paper) [W*m^-2*K^-1] 
k_dLv = k_dEff/l_v; % calibration constant (see sect 4.1 of dihiraj's DEID paper) [kg*m^-2*K^-1*s^-1] 
c_melt = 6.69e-05; % melting calibration constant (see sect. 4.1 of dhiraj density paper) [m/s*K]  

% Eqn. (13) in Dhiraj's density paper: c = (L_vv) / (L_ff*C_melt)
% c = hf_rho_coeff [K*s/m]
hf_rho_coeff_calc = (l_v)/ (l_f*c_melt); 
hf_rho_coeff = 1.01e05; % this is the value pulled from Dhiraj's paper using l_vv and l_ff  


%% Move to working directory and identify most recent video file
% Sets path for .m to run in:
cd(working_dir) 
% Get a list of all video filesin this directory
directory = dir("*.avi");
% Find the most recent .avi file in this directory 
[~,idx] = max([directory.datenum]);
latest_file =  directory(idx).name;
%% When testing: 
% % Move to working directory and identify video files 
% % Sets path for .m to run in:
% cd(working_dir) 
% % Get a list of all files and folders in this directory
% directory = dir(".");
% % Initialize an empty cell array to store the file names
% file_names = {};
% % Loop through each item in the input directory
% for file_i = 1:length(directory)
%     % Check if the item is a file (not a folder) and if it ends with .avi
%     if ~directory(file_i).isdir && endsWith(directory(file_i).name, '.avi', 'IgnoreCase', true)
%         % Get the name of the file and append it to the list
%         file_names{end+1} = directory(file_i).name;
%     end
% end

%% Initialize Output Tables
% Frame by Frame output table 
% fbf_col_names = {'Time', 'SWE_mm'};
% fbf_col_types = {'datetime', 'double'};
% fbf_output_table = table('Size', [0, length(fbf_col_names)], ...
%                          'VariableNames', fbf_col_names, ...
%                          'VariableTypes', fbf_col_types);
% Time series output table
ts_col_names = {'Time', 'Complexity', 'SDI', 'Mass', 'Volume', 'Diameter', 'Surface Area', 'Void Space', 'Density [kg/m^3]', 'SWE [mm]','Snow [mm]'};
ts_col_types = {'datetime', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
ts_output_table = table('Size', [0, length(ts_col_names)], ...
                         'VariableNames', ts_col_names, ...
                         'VariableTypes', ts_col_types);
% Particles output table
particle_col_names = {'Time', 'Terminal_Velocity', 'Complexity', 'SDI', 'Mass', 'Volume HFD', 'Volume Sphere', 'Density HFD', 'Density Sphere', 'Diameter', 'Surface Area', 'Void Space', 'PBP SWE (mm)', 'FBF SWE (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'SWE Factor', 'Missing Data'};
particle_col_types = {'datetime', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double' ,'double', 'double', 'double', 'double', 'double', 'logical'};
particle_output_table = table('Size', [0, length(particle_col_names)], ...
                         'VariableNames', particle_col_names, ...
                         'VariableTypes', particle_col_types);
% DEID summary table 
summary_col_names = {'Time', 'Complexity', 'SDI', 'HFD Density (kg*m^-3)', 'PBP SWE (mm)', 'FBF SWE (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'Hot Plate Area', 'SWE Factor'};
summary_col_types = {'datetime', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
DEID_summary_table = table('Size', [0, length(summary_col_names)], ...
                         'VariableNames', summary_col_names, ...
                         'VariableTypes', summary_col_types);
% % Diagnostic output table
% diag_col_names = {'Filename', 'Start_Time','End_Time', 'SWE_Factor', 'Num_Particles'};
% diag_col_types = {'string', 'datetime', 'datetime', 'double', 'double'};
% diag_output_table = table('Size', [0, length(diag_col_names)], ...
%                          'VariableNames', diag_col_names, ...
%                          'VariableTypes', diag_col_types);
% diag_table = table('Size', [1, length(diag_col_names)], ...
%                          'VariableNames', diag_col_names, ...
%                          'VariableTypes', diag_col_types);

%% Begin DEID video processing:

% for file_i = 1:length(file_names)

    filename = latest_file;
    % filename = file_names{file_i};
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
    % diag_table.Filename = filename;
    % diag_table.Start_Time = vid_start_time;
    % diag_table.End_Time = vid_end_time;

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
    hp_area = size(frame_gray_cropped,1) * size(frame_gray_cropped,2) * pix_to_m2_conversion;    
    % hp_area(file_i) = hp_area_temp; % Hotplate Area [m^2]
    h_mass_fbf = (k_dLv*sum_h_area_times_dt) / vid_fps; % mass evaporates in each frame
    SWE_FBF_mm = h_mass_fbf / hp_area;
    % Find the minimum SWE in all frames within a video, and subtract from
    % SWE (way of handling residue) 
    fbf_SWE_min = min(SWE_FBF_mm);
    SWE_FBF_mm = SWE_FBF_mm - fbf_SWE_min;
    SWE_fbf_accumulation = sum(SWE_FBF_mm);
    % SWE_fbf_accumulation_noSub = sum(SWE_FBF_mm)
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
    h_energy_per_time = h_mass_pbp * l_vv ./ (h_max_area .* h_evap_time); % Heat flux method: energy per unit area per time
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
    % for current testing, ben is abbreivating the list of variables
    % stored:

    pbp_table_particles = table(h_initial_time', h_evap_time', h_mass_pbp', h_diam',  ...
        h_max_area', h_max_circ_area', h_rho_sph', h_rho_hfd',  ... 
        h_vol_hfd', h_sph_vol');
    pbp_table_particles.Properties.VariableNames = {'initial_time', 'evap_time', 'mass', 'diam',  ...
        'max_area', 'max_circ_area', 'rho_sph', 'rho_hfd', ... 
        'vol_hfd', 'vol_sph'};

    % % Handles PBP data
    % pbp_table_particles = table(h_initial_time', h_evap_time', h_mass_pbp', h_diam',  ...
    %     h_max_area', h_max_circ_area', h_rho_sph', h_rho_hfd',  ... 
    %     h_vol_hfd', h_water_eq_diameter', h_max_maj_axis', h_max_min_axis', ...
    %     h_energy_per_time', h_delta_temp_max', h_delta_temp_mean');
    % pbp_table_particles.Properties.VariableNames = {'initial_time', 'evap_time', 'mass', 'diam',  ...
    %     'max_area', 'max_circ_area', 'rho_sph', 'rho_hfd', ... 
    %     'vol_hfd', 'water_eq_diameter', 'max_maj_axis', 'max_min_axis', ...
    %     'energy_per_time', 'delta_temp_max', 'delta_temp_mean'};
    
    pbp_table_particles = table2timetable(pbp_table_particles); % Convert to timetable 
    pbp_table_particles = sortrows(pbp_table_particles, 'initial_time'); % Sort by time 
    % diag_table.Num_Particles = size(pbp_table_particles,1);

    % Filters data to find where 0 < mass < .005 to omit residue on plate
    g1 = find((pbp_table_particles.mass > 0 &... 
        pbp_table_particles.mass < residue_filter) &...
        (pbp_table_particles.evap_time > evapTime_min &... 
        pbp_table_particles.evap_time < evapTime_max)); 
        
    pbp_table_particles = pbp_table_particles(g1,:);
    %% Post Processing Starts Here! 
    % Processes PBP data if particles were present in video file
    if height(pbp_table_particles)>100

        % Calculate terminal velocity using terminalVelocity function
        pbp_table_particles.terminal_vel = terminalVelocity(pbp_table_particles.max_area, pbp_table_particles.max_circ_area, pbp_table_particles.mass);
        
        % Calculate 'void space'
        pbp_table_particles.void_space = (pbp_table_particles.max_area - pbp_table_particles.max_circ_area).^(1/2); 
        
        for i = 1:length(pbp_table_particles.void_space)
            if pbp_table_particles.void_space(i) > 0
                pbp_table_particles.void_space(i) = pbp_table_particles.void_space(i);
            else
                pbp_table_particles.void_space(i) = 0;
            end
        end

        % Calculate SDI, Cx  
        pbp_table_particles.surface_area_eq = ((9*pi)/16)^(1/3) * (pbp_table_particles.mass ./ rho_water).^(2/3); % Surface area of equivalent water droplet (See POF Singh et al. 2023)
        pbp_table_particles.complexity = pbp_table_particles.max_circ_area ./ pbp_table_particles.max_area; % 'Complexity' (See CRST Morrison et al. 2023)
        pbp_table_particles.sdi = pbp_table_particles.max_area ./ pbp_table_particles.surface_area_eq ; % 'SDI' (See CRST Morrison et al. 2023)

        % Resamples time series at desired interval
        fbf_table_retimed = retime(fbf_table_raw, 'regular', 'sum', 'TimeStep', time_step);
        fbf_table_retimed.SWE_FBF_accum_mm = cumsum(fbf_table_retimed.SWE_FBF_mm);

        %% Total SWE for all PBP data
        pbp_table_particles.SWE_PBP_mm = c1 * pbp_table_particles.mass ./ (rho_water * hp_area); % [mm]
        pbp_table_particles.SWE_PBP_accum_mm = cumsum(pbp_table_particles.SWE_PBP_mm);  % [mm]
     
        % Finds difference factor between FBF SWE and PBP SWE and adjusts PBP SWE
        swe_factor = fbf_table_retimed.SWE_FBF_accum_mm(end) / pbp_table_particles.SWE_PBP_accum_mm(end);
        % if swe factor is abnormally high, adjust it to account for
        % residue:
        if swe_factor > 1.95
            swe_factor = 1.95; 
        end

        % Adjusts PBP SWE
        pbp_table_particles.SWE_FBF_mm = pbp_table_particles.SWE_PBP_mm * swe_factor; 
        pbp_table_particles.SWE_FBF_accum_mm = cumsum(pbp_table_particles.SWE_FBF_mm);
        
        % Calculate SWE factor for all PBP data:
        pbp_table_particles.SWEfactor = pbp_table_particles.SWE_FBF_mm ./ pbp_table_particles.SWE_PBP_mm; 
        %% Total Snow for all PBP data
        pbp_table_particles.snow_PBP_mm = rho_water * (pbp_table_particles.SWE_PBP_mm ./ pbp_table_particles.rho_hfd); % [mm]
        pbp_table_particles.snow_PBP_acc_mm = cumsum(pbp_table_particles.snow_PBP_mm); % [mm]
        pbp_table_particles.snow_FBF_mm = rho_water * (pbp_table_particles.SWE_FBF_mm ./ pbp_table_particles.rho_hfd); % [mm]
        pbp_table_particles.snow_FBF_acc_mm = cumsum(pbp_table_particles.snow_FBF_mm); % [mm]
        
        %% Appends PARTICLE data for single video to output table
        % selects a subset of output variables to be exported:
        particle_output_table = pbp_table_particles(:, {'terminal_vel', 'complexity', 'sdi', 'mass', 'vol_hfd', 'vol_sph', 'rho_hfd', 'rho_sph', 'diam', 'surface_area_eq', 'void_space', 'SWE_PBP_mm', 'SWE_FBF_mm', 'snow_PBP_mm', 'snow_FBF_mm', 'SWEfactor'});
        particle_output_table = timetable2table(particle_output_table);
        % add a flag column to distinguish missing .avi data:
        particle_output_table.missing_data = false(height(particle_output_table),1); 
        particle_output_table.Properties.VariableNames = particle_col_names;
        
        % Loop through each variable and replace NaNs with zeros.
        % NaNs are showing up for some videos because the pbp and fbf tables
        % dont have the same number of time stamps. Replacing Nans is a quick
        % fix but should be further explored.
        % particle by particle
        output_names = particle_output_table.Properties.VariableNames;
        for i = 2:length(output_names)
            particle_output_table.(output_names{i})(isnan(particle_output_table.(output_names{i}))) = 0;
        end

         % Appends output
         % particle_output_table_all = [particle_output_table_all; particle_output_table_vid(:, particle_col_names)];

        %% Find an average/sum of previous two minutes of data to fill time gap between .avi files 

        % Go back two minutes and capture corresponding data:
        final_time = particle_output_table.Time(end);
        prev_time = final_time - minutes(2);
        prev_data = particle_output_table(particle_output_table.Time >= prev_time, :);

        % take the mean or sum of each value corresponding to the captured data:
        termVelocityRow = mean(prev_data.Terminal_Velocity);
        cxRow = mean(prev_data.Complexity);
        sdiRow = mean(prev_data.SDI);
        massRow = mean(prev_data.Mass);
        volumeHFDrow = mean(prev_data.('Volume HFD'));
        volumeSPHrow = mean(prev_data.('Volume Sphere'));
        densityHFDrow = massRow / volumeHFDrow;
        densitySPHrow = massRow / volumeSPHrow;
        diamRow = mean(prev_data.Diameter);
        surAreaRow = mean(prev_data.('Surface Area'));
        voidSpaceRow = mean(prev_data.('Void Space'));
        PBPsweRow = 1000 * massRow/ (rho_water * hp_area);
        FBFsweRow = PBPsweRow*mean(prev_data.('SWE Factor'));
        PBPsnowRow = rho_water * (PBPsweRow ./ densityHFDrow);
        FBFsnowRow = rho_water * (FBFsweRow ./ densityHFDrow);
        SWEfactorRow = mean(prev_data.('SWE Factor')); 
        % compile into a table
        new_row = array2table([termVelocityRow, cxRow, sdiRow, massRow, ...
            volumeHFDrow, volumeSPHrow, densityHFDrow, densitySPHrow, diamRow, surAreaRow, voidSpaceRow,...
            PBPsweRow, FBFsweRow PBPsnowRow, FBFsnowRow, SWEfactorRow], 'VariableNames', particle_output_table.Properties.VariableNames(2:(end-1)));
        % assign a time and logical value for missing data to the new row:
        new_row.('Missing Data') = true;
        new_row.Time = final_time + seconds(5);

        % now append new row to particle_output_table:
        particle_output_table = [particle_output_table; new_row];

        %% Averaging technique starts here! 
        % particle_output = timetable2table(particle_output_table);
        particle_output_table = table2timetable(particle_output_table)  
          
        %% Computes averaged and summed PBP data using retime function

        % Averages PBP data 
        avg_cols = {'Complexity', 'SDI', 'Diameter', 'Surface Area', 'Void Space', 'SWE Factor'};
        avg_table = retime(particle_output_table(:, avg_cols), 'regular', 'mean', 'TimeStep', time_step);
        % Sums PBP data
        sum_cols = {'Mass', 'Volume HFD'};
        sum_table = retime(particle_output_table(:, sum_cols), 'regular', 'sum', 'TimeStep', time_step);
        % Join tables
        pbp_table_retimed = horzcat(avg_table, sum_table);
        
        % calculate density:
        pbp_table_retimed.density = pbp_table_retimed.Mass ./ pbp_table_retimed.('Volume HFD');   % HFD Density method [kg/m^3]
        % calculate SWE:
        pbp_table_retimed.SWE_PBP_mm = 1000 * pbp_table_retimed.Mass ./ (rho_water * hp_area(1)); % [mm]
        pbp_table_retimed.SWE_FBF_mm = pbp_table_retimed.SWE_PBP_mm.*pbp_table_retimed.('SWE Factor'); % [mm]
        % pbp_table_retimed.SWE_PBP_accum_mm = cumsum(pbp_table_retimed.SWE_PBP_mm);  % [mm]
        pbp_table_retimed.SWE_FBF_accum_mm = cumsum(pbp_table_retimed.SWE_FBF_mm);  % [mm]
        % calculate snow:
        pbp_table_retimed.Snow_mm = rho_water * pbp_table_retimed.SWE_FBF_mm./ pbp_table_retimed.density; % [mm]
        pbp_table_retimed.Snow_accum_mm = cumsum(pbp_table_retimed.Snow_mm); % [mm]
        
        %% Appends AVERAGED TIME SERIES data for single video to output table
        % video_output_table = synchronize(fbf_table_retimed, pbp_table_retimed);
        % Selects a subset of output variables to be exported
        ts_output_table = pbp_table_retimed(:, {'Complexity', 'SDI', 'Mass', 'Volume HFD', 'Diameter', 'Surface Area', 'Void Space', 'density', 'SWE_FBF_mm','Snow_mm'});
        ts_output_table = timetable2table(ts_output_table);
        ts_output_table.Properties.VariableNames = ts_col_names;
        
        % Loop through each variable and replace NaNs with zeros.
        % NaNs are showing up for some videos because the pbp and fbf tables
        % dont have the same number of time stamps. Replacing Nans is a quick
        % fix but should be further explored.
        output_names = ts_output_table.Properties.VariableNames;
        for i = 2:length(output_names)
            ts_output_table.(output_names{i})(isnan(ts_output_table.(output_names{i}))) = 0;
        end
        %% Sorts table by time and handles duplicates 
        % particle_output_table = table2timetable(sortrows(particle_output_table, 'Time'));
        ts_output_table = sortrows(ts_output_table, 'Time');
        ts_output_table = table2timetable(ts_output_table); 
        
        %% Cumulatively sums data for SWE and Snow totals
        % Averaged data 
        ts_output_table.('SWE Total [mm]') = cumsum(ts_output_table.('SWE [mm]'));
        ts_output_table.('Snow Total [mm]') = cumsum(ts_output_table.('Snow [mm]'));
        ts_output_table.('Snow Total [in]') = ts_output_table.('Snow Total [mm]') * mm_to_inches;

        %% Now create a summary table with just total SWE, Snow, and average density per .avi 
        
        DEID_summary_table = timetable(time_series(end));
        DEID_summary_table.cx = mean(particle_output_table.Complexity);
        DEID_summary_table.sdi = mean(particle_output_table.SDI);
        DEID_summary_table.rho = sum(particle_output_table.Mass) / sum(particle_output_table.('Volume HFD'));
        DEID_summary_table.pbpSWE = sum(particle_output_table.('PBP SWE (mm)'));
        DEID_summary_table.fbfSWE = sum(particle_output_table.('FBF SWE (mm)'));
        DEID_summary_table.pbpSnow = sum(particle_output_table.('PBP Snow (mm)'));
        DEID_summary_table.fbfSnow = sum(particle_output_table.('FBF Snow (mm)'));
        DEID_summary_table.hotPlateArea = hp_area; 
        DEID_summary_table.SWEfactor = swe_factor; 
        DEID_summary_table = timetable2table(DEID_summary_table);
        DEID_summary_table.Properties.VariableNames = summary_col_names; 
        DEID_summary_table = table2timetable(DEID_summary_table);
        
        %% Save processed output data for all video files present
        
        % Get folder name and saves output as 'folder name'.csv:
        start_time = datestr(ts_output_table.Time(1), 'yyyy-mm-dd_HH-MM-ss');
        current_time = datestr(now, 'mm-dd-yyyy');
        
        % Writes out all particle data table:
        writetimetable(particle_output_table, [output_dir, '\DEID_Particle_', start_time, '.csv'], 'Delimiter', ',');
        
        % Writes out time averaged data table: 
        writetimetable(ts_output_table, [output_dir, '\DEID_TS_', start_time, '.csv'], 'Delimiter', ',');
        
        % Writes out DEID summary table for current storm:
        writetimetable(DEID_summary_table, [output_dir, '\DEID_totals.csv'], 'WriteMode', 'append', 'Delimiter', ',');
        
        % Writes out DEID summary table for each day of storm:
        writetimetable(DEID_summary_table, [output_dir, '\DEID_totals_', current_time, '.csv'], 'WriteMode', 'append', 'Delimiter', ',');

        % Converts DEID summary table to json:
        % jsonPath = 'D:\Atwater\DEID_totals.json';
        % jsonTable = jsonencode(timetable2table(DEID_summary_table));
        % fid = fopen(jsonPath, 'w')
        % if fid == -1
        %     error('cannot open file for writing');
        % end
        % fwrite(fid, jsonTable, 'char');
        % fclose(fid);
        
        %% Send each .csv file to chpc
        
        % DEID storm summary table:
        % Construct filename:
        summary_file = sprintf('%s\\DEID_totals.csv', output_dir);
        % Construct the SCP command
        scpCommand_summary = sprintf('"C:\\Program Files\\PuTTY\\pscp.exe" -pw %s "%s" %s@%s:%s', '"Sc0tchT@p3!"', summary_file, 'u6022893', 'notchpeak2.chpc.utah.edu', '/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025');
        status_summary = system(scpCommand_summary);

        % Check if the command was successful
        if status_summary == 0
            disp('DEID_totals.csv transferred succesfully');
        else
            disp('Error transferring DEID_totals.csv');
        end

        % DEID storm summary json:
        % Construct the SCP command
        % scpCommand_json = sprintf('"C:\\Program Files\\PuTTY\\pscp.exe" -pw %s "%s\DEID_totals.json" %s@%s:%s', '"Sc0tchT@p3!"', output_dir, 'u6022893', 'notchpeak2.chpc.utah.edu', '/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025');
        % status_json = system(scpCommand_json);

        % Check if the command was successful
        % if status_json == 0
        %     disp('DEID_totals.json transferred succesfully');
        % else
        %     disp('Error transferring DEID_json.csv');
        % end

        % DEID particle file:
        % Construct filename:
        particle_file = sprintf('%s\\DEID_Particle_%s.csv', output_dir, start_time);
        % Construct the SCP command
        scpCommand_particle = sprintf('"C:\\Program Files\\PuTTY\\pscp.exe" -pw %s "%s" %s@%s:%s', '"Sc0tchT@p3!"', particle_file, 'u6022893', 'notchpeak2.chpc.utah.edu', '/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025');
        status_particle = system(scpCommand_particle);

        % Check if the command was successful
        if status_particle == 0
            disp(['File ', particle_file, ' transferred succesfully.']);
        else
            disp(['Error transferring ', particle_file, '.']);
        end

        % DEID time series file:
        % Construct filename:
        TS_file = sprintf('%s\\DEID_TS_%s.csv', output_dir, start_time);
        % Construct the SCP command
        scpCommand_TS = sprintf('"C:\\Program Files\\PuTTY\\pscp.exe" -pw %s "%s" %s@%s:%s', '"Sc0tchT@p3!"', TS_file, 'u6022893', 'notchpeak2.chpc.utah.edu', '/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025');
        % Execute the SCP command using the system function
        % system('cd C:\Program Files\PuTTY\')
        status_TS = system(scpCommand_TS);

        % Check if the command was successful
        if status_TS == 0
            disp(['File ', TS_file, ' transferred succesfully.']);
        else
            disp(['Error transferring ', TS_file, '.']);
        end

        % DEID day summary file:
        % Construct filename:
        daySummary_file = sprintf('%s\\DEID_totals_%s.csv', output_dir, current_time);
        % Construct the SCP command
        scpCommand_daySummary = sprintf('"C:\\Program Files\\PuTTY\\pscp.exe" -pw %s "%s" %s@%s:%s', '"Sc0tchT@p3!"', daySummary_file, 'u6022893', 'notchpeak2.chpc.utah.edu', '/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025');
        status_daySummary = system(scpCommand_daySummary);
         
        % Check if the command was successful
        if status_daySummary == 0
             disp(['File ', daySummary_file, ' transferred succesfully.']);
        else
             disp(['Error transferring ', daySummary_file, '.']);
        end

    else
        
        %% Create blank tables:
        % summary table:
        DEID_summary_table = timetable(time_series(end));
        DEID_summary_table.cx = 0;
        DEID_summary_table.sdi = 0;
        DEID_summary_table.rho = 0;
        DEID_summary_table.pbpSWE = 0;
        DEID_summary_table.fbfSWE = 0;
        DEID_summary_table.pbpSnow = 0;
        DEID_summary_table.fbfSnow = 0;
        DEID_summary_table.hotPlateArea = hp_area; 
        DEID_summary_table.SWEfactor = swe_factor; 
        DEID_summary_table = timetable2table(DEID_summary_table);
        DEID_summary_table.Properties.VariableNames = summary_col_names; 
        DEID_summary_table = table2timetable(DEID_summary_table);
        % particle table:
        particle_output_table = table('Size', [0, length(particle_col_names)], ...
                         'VariableNames', particle_col_names, ...
                         'VariableTypes', particle_col_types);
        % time series table:
        ts_output_table = table('Size', [0, length(ts_col_names)], ...
                         'VariableNames', ts_col_names, ...
                         'VariableTypes', ts_col_types);
        % Convert all tables to time tables if they are tables:
        if istable(particle_output_table)
            particle_output_table = table2timetable(particle_output_table);
            ts_output_table = table2timetable(ts_output_table);
        end
        % Get folder name and saves output as 'folder name'.csv:
        start_time = datestr(time_series(1), 'yyyy-mm-dd_HH-MM-ss');
        current_time = datestr(now, 'mm-dd-yyyy');
        
        % Writes out all particle data table:
        writetimetable(particle_output_table, [output_dir, '\DEID_Particle_', start_time,'.csv'], 'Delimiter', ',');
        
        % Writes out time averaged data table: 
        writetimetable(ts_output_table, [output_dir, '\DEID_TS_', start_time, '.csv'], 'Delimiter', ',');
        
        % Writes out DEID summary table for current storm:
        writetimetable(DEID_summary_table, [output_dir, '\DEID_totals.csv'], 'WriteMode', 'append', 'Delimiter', ',');
        
        % Writes out DEID summary table for each day of storm:
        writetimetable(DEID_summary_table, [output_dir, '\DEID_totals_', current_time, '.csv'], 'WriteMode', 'append', 'Delimiter', ',');

        % Converts DEID summary table to json:
        % jsonPath = 'D:\Atwater\DEID_totals.json';
        % jsonTable = jsonencode(timetable2table(DEID_summary_table));
        % fid = fopen(jsonPath, 'w')
        % if fid == -1
        %     error('cannot open file for writing');
        % end
        % fwrite(fid, jsonTable, 'char');
        % fclose(fid);
        
        %% Send each .csv file to chpc
        
        % DEID storm summary table:
        % Construct filename:
        summary_file = sprintf('%s\\DEID_totals.csv', output_dir);
        % Construct the SCP command
        scpCommand_summary = sprintf('"C:\\Program Files\\PuTTY\\pscp.exe" -pw %s "%s" %s@%s:%s', '"Sc0tchT@p3!"', summary_file, 'u6022893', 'notchpeak2.chpc.utah.edu', '/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025');
        status_summary = system(scpCommand_summary);

        % Check if the command was successful
        if status_summary == 0
            disp('DEID_totals.csv transferred succesfully');
        else
            disp('Error transferring DEID_totals.csv');
        end

        % DEID storm summary json:
        % Construct the SCP command
        % scpCommand_json = sprintf('"C:\\Program Files\\PuTTY\\pscp.exe" -pw %s "%s\DEID_totals.json" %s@%s:%s', '"Sc0tchT@p3!"', output_dir, 'u6022893', 'notchpeak2.chpc.utah.edu', '/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025');
        % status_json = system(scpCommand_json);

        % Check if the command was successful
        % if status_json == 0
        %     disp('DEID_totals.json transferred succesfully');
        % else
        %     disp('Error transferring DEID_json.csv');
        % end

        % DEID particle file:
        % Construct filename:
        particle_file = sprintf('%s\\DEID_Particle_%s.csv', output_dir, start_time);
        % Construct the SCP command
        scpCommand_particle = sprintf('"C:\\Program Files\\PuTTY\\pscp.exe" -pw %s "%s" %s@%s:%s', '"Sc0tchT@p3!"', particle_file, 'u6022893', 'notchpeak2.chpc.utah.edu', '/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025');
        status_particle = system(scpCommand_particle);

        % Check if the command was successful
        if status_particle == 0
            disp(['File ', particle_file, ' transferred succesfully.']);
        else
            disp(['Error transferring ', particle_file, '.']);
        end

        % DEID time series file:
        % Construct filename:
        TS_file = sprintf('%s\\DEID_TS_%s.csv', output_dir, start_time);
        % Construct the SCP command
        scpCommand_TS = sprintf('"C:\\Program Files\\PuTTY\\pscp.exe" -pw %s "%s" %s@%s:%s', '"Sc0tchT@p3!"', TS_file, 'u6022893', 'notchpeak2.chpc.utah.edu', '/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025');
        % Execute the SCP command using the system function
        % system('cd C:\Program Files\PuTTY\')
        status_TS = system(scpCommand_TS);

        % Check if the command was successful
        if status_TS == 0
            disp(['File ', TS_file, ' transferred succesfully.']);
        else
            disp(['Error transferring ', TS_file, '.']);
        end

        % DEID day summary file:
        % Construct filename:
        daySummary_file = sprintf('%s\\DEID_totals_%s.csv', output_dir, current_time);
        % Construct the SCP command
        scpCommand_daySummary = sprintf('"C:\\Program Files\\PuTTY\\pscp.exe" -pw %s "%s" %s@%s:%s', '"Sc0tchT@p3!"', daySummary_file, 'u6022893', 'notchpeak2.chpc.utah.edu', '/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025');
        status_daySummary = system(scpCommand_daySummary);
         
        % Check if the command was successful
        if status_daySummary == 0
             disp(['File ', daySummary_file, ' transferred succesfully.']);
        else
             disp(['Error transferring ', daySummary_file, '.']);
        end
    end

[~, parent_dir, ~] = fileparts(pwd);
disp(['Saved Output for: ', parent_dir])

% end 

exit 
