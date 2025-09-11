%% DEID AVI File Processing Code
% Outputs particle by particle csclear,v file to be post processed by
% DEID_AVI_Processor.m
% AUTHOR : Dhiraj Singh, Benjamin Silberman, Travis Morrison, Alex Blackmer
                                       
clear, clc, close all
%% set filepath, output directory, and file name for saving  
working_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/DEC/DEC_all';
output_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/test';
storm_output = '_DEC_2023';

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
% DEID specific parameters:
residue_filter = 0.005; % max weight of snowflake to process [kg]
evapTime_min = 1/15; % minimum time a snowflake has to appear on hotplate to be processed
evapTime_max = 10; % maximum time a snowflake can appear on hotplate to be processed
colorbar_image_indexes = [1 1 384 288]; % location of colorbar in pixel locations
crop_index = 43; % use this to specify indices to crop out kapton tape
colorbar_kapton_image_indexes = [1 (colorbar_image_indexes(2)+crop_index) 383 (colorbar_image_indexes(4)-crop_index)]; % location of Kapton tape in pixel locations
colorbar_max_temp = 145; % max temperature set in colorbar on the physical screen of the tir software
min_thres = 70; % minimum threshold number in image accepted rbg ([0 255]) scale
sort_threshold = 20; % this is the RMS threshold between succesive images of snowflakes used in the sortPostitions_v2
min_h_size = 10; % minumum hydrometeor size in pixels 
minimum_drop_life = 0; % minimum number of frames a drop has to be visable to be processed
maximum_drop_life = 150; % maximum number of frames a drop has to be visable to be processed
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
file_names = file_names(1:5);

% %% loop through each file_name to get start and end times:
% % this is not totally necessary and can be quite time extensive, 
% % but it is helpful when testing
% vid_end_time = NaT(length(file_names), 1);
% vid_start_time = NaT(length(file_names), 1);
% 
% parfor file_i = 1:length(file_names)
%     filename = file_names{file_i};
%     vid=VideoReader(filename);
%     % get metaData for video processing
%     vid_dir = dir(filename);
%     vid_fps = vid.FrameRate;
%     num_frames = vid.NumFrames;
%     % create a time series of date times starting from (video end time - video duration) and ending at the video end time
%     vid_end_time(file_i) = datetime([vid_dir.date]);
%     time_series = datetime(vid_end_time(file_i) - (0:num_frames) * seconds(1/vid_fps), 'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
%     time_series = flip(time_series);  % flips time series to be chronologically ordered
%     vid_start_time(file_i) = datetime(time_series(1));
% end
% 
% start_end_time_table = table(file_names', vid_start_time, vid_end_time); 
% start_end_time_table.length = start_end_time_table.vid_end_time - start_end_time_table.vid_start_time; 
% start_end_time_table.Properties.VariableNames{1} = 'file_name'; 

%% initialize Output Tables
% time series output table
ts_col_names = {'Time', 'Complexity', 'SDI', 'Mass', 'Volume', 'Diameter', 'Surface Area', 'Void Space', 'Density [kg/m^3]', 'SWE [mm]','Snow [mm]'};
ts_col_types = {'datetime', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
ts_output_table = table('Size', [0, length(ts_col_names)], ...
                         'VariableNames', ts_col_names, ...
                         'VariableTypes', ts_col_types);
% particles output table
particle_col_names = {'Time', 'Terminal_Velocity', 'Complexity', 'SDI', 'Evap Time', 'Temp Diff', 'Snowflake Height', 'Mass', 'Volume HFD', 'Volume Sphere', 'Density HFD', 'Density Sphere', 'Diameter', 'Surface Area', 'Void Space', 'PBP SWE (mm)', 'FBF SWE (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'SWE Factor', 'Missing Data'};
particle_col_types = {'datetime', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double' ,'double', 'double', 'double', 'double', 'double', 'logical'};
particle_output_table_all = table('Size', [0, length(particle_col_names)], ...
                         'VariableNames', particle_col_names, ...
                         'VariableTypes', particle_col_types);
% DEID summary table 
summary_col_names = {'Time', 'Complexity', 'SDI', 'HFD Density (kg*m^-3)', 'PBP SWE (mm)', 'FBF SWE (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'Hot Plate Area', 'SWE Factor'};
summary_col_types = {'datetime', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
DEID_summary_table = table('Size', [0, length(summary_col_names)], ...
                         'VariableNames', summary_col_names, ...
                         'VariableTypes', summary_col_types);

%% begin DEID video processing:

% Parfor loop parallelizes processing by distributing each video file to a
% Matlab worker on each CPU core. 

% preallocate variables:
pbp_table_particles_cell = cell(length(file_names),1);
vid_end_time = NaT(length(file_names), 1);
vid_start_time = NaT(length(file_names), 1);
vid_end_start_diff = NaN(length(file_names)-1, 1);
final_particle_output_table = table(); 
hp_area = NaN(1, length(file_names)); 
SWE_factor = NaN(1, length(file_names)); 
fbf_SWE_min = NaN(1, length(file_names)); 

 for file_i = 1:length(file_names)
    filename = file_names{file_i};
    disp(['Processing File: ', filename])
    vid=VideoReader(filename);
    
    % get metaData for video processing 
    vid_dir = dir(filename);
    vid_fps = vid.FrameRate;
    num_frames = vid.NumFrames;

    % gets video start and end date to construct time series
    vid_end_time(file_i) = datetime([vid_dir.date]);
    
    % create a time series of date times starting from (video end time - video duration) and ending at the video end time
    time_series = datetime(vid_end_time(file_i) - (0:num_frames) * seconds(1/vid_fps), 'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
    time_series = flip(time_series);  % Flips time series to be chronologically ordered
    vid_start_time(file_i) = datetime(time_series(1));

    % preallocate cell for hydrometor data
    h_data = cell(num_frames,1); 

    %% "frame by frame method"; this is how we obtain SWE for each .avi file 
    % preallocate variables saved in loop for speed:
    plate_temp = nan(num_frames,1);
    sum_h_area_times_dt = nan(num_frames,1);
    % enter loop to process images: 
    for frame_ii = 1:num_frames     
        frame = read(vid, frame_ii);  
        frame_gray = im2gray(frame); % convert frame of interest to gray scale
        frame_gray_cropped_wKapton = imcrop(frame_gray, colorbar_image_indexes);% crop out colorbar
        plate_temp(frame_ii) = max(max(double(frame_gray_cropped_wKapton))); % this assumes max temperature in image is the plate temperature with Kapton tape 
        frame_gray_cropped = imcrop(frame_gray, colorbar_kapton_image_indexes); % back to orginal grayscale image... now remove colorbar and kapton tape from image
        frame_filtered = frame_gray_cropped > min_thres; % removed below min threshold, on rbg ([0, 255]) scale 
        frame_filtered_filled = imfill(frame_filtered, 'Holes'); % clean up Hydrometeors
        % Get hydrometeor properties: 
        h_geo_prop = regionprops(frame_filtered_filled, 'MajorAxisLength', 'MinorAxisLength', 'Centroid', 'Area','BoundingBox'); % returns the centroid, the area of each blob, and the bounding box (left, top, width, height).
        % if no properties are found, go to next frame: 
        if (isempty(h_geo_prop))
            continue;
        end
        % build Hydrometeor property matrices from regionprops values: 
        h_bounding_box = cat(1,h_geo_prop.BoundingBox); % concat all values to build matrix of bounding box indices
        h_centroid = round(cat(1, h_geo_prop.Centroid)); % concat all values to build matrix of centroids
        h_major = cat(1,h_geo_prop.MajorAxisLength);
        h_minor = cat(1,h_geo_prop.MinorAxisLength);
        h_area_pix = cat(1, h_geo_prop.Area); % concat all values to build matrix of Hydrometeor areas in pixels
        % convert Hydrometeor areas to m^2: 
        h_area = h_area_pix .* pix_to_m2_conversion; 
        h_major_axis = h_bounding_box(:, 3) * sqrt(pix_to_m2_conversion); % hydrometeor major axis
        h_minor_axis = h_bounding_box(:, 4) * sqrt(pix_to_m2_conversion); % hydrometeor minor axis
        h_major = h_major * pix_to_m_conversion;
        h_minor = h_minor * pix_to_m_conversion;
        % get the Hydrometeor ellipse area:
        h_elipse_area = h_major_axis .* h_minor_axis;
        % difference in temperature of each centroid and the plate:
        h_centroid_i = sub2ind(size(frame_gray_cropped), h_centroid(:, 2), h_centroid(:, 1)); % find the index of the centriods in orginal image
        h_centroid_values = double(frame_gray_cropped(h_centroid_i)); % intensities of centroid pixels of snow
        plate_h_dtemp = colorbar_max_temp - (h_centroid_values .* (colorbar_max_temp / plate_temp(frame_ii)));
        % product of hydrometeor area with the temp difference:
        h_area_times_dtemp = h_area .* plate_h_dtemp;         
        % sum the product of individual area and temp. diff in each frame:
        % **this is how we obtain hydrometeor mass using fbf method**
        sum_h_area_times_dt(frame_ii)=sum(h_area_times_dtemp);         
        % Build large matrix of Hydrometeor data
        h_data{frame_ii} = cat(2, h_centroid, h_area, plate_h_dtemp, h_elipse_area, h_major_axis, h_minor_axis); 
    end

    % frame by frame SWE calculation
    hp_area(file_i) = size(frame_gray_cropped,1) * size(frame_gray_cropped,2) * pix_to_m2_conversion; % hotplate area     
    h_mass_fbf = (k_dLv*sum_h_area_times_dt) / vid_fps; % total mass evaporates in each frame
    SWE_fbf = h_mass_fbf / hp_area(file_i);
    SWE_fbf = SWE_fbf - min(SWE_fbf); 
    % find the minimum SWE in all frames within a video, and subtract from
    % SWE (way of handling residue) 
    fbf_SWE_min(file_i) = min(SWE_fbf(SWE_fbf ~=0));
    SWE_fbf = SWE_fbf - fbf_SWE_min(file_i); 
    time_series_fbf = time_series(1:length(SWE_fbf));
    
    %% sorting one frame to others frame ~ data cleaning of some sorts: 
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
    
    %% isolating the variables and put them into a matrix to work with:
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
    % preallocate variables:
    h_delta_time = {}; % hydrometeor time to melt and evaporate
    h_delta_temp_max = {}; % hydrometeor centroid's max intensity
    h_delta_temp_mean = {}; % hydrometeor centroid's mean intensity
    h_max_area = {}; % hydrometeor maximum area
    h_max_circ_area = {};% hydrometeor Max circumscribed area
    h_max_maj_axis = {};% hydrometeor Max major axis
    h_max_min_axis = {};% hydrometeor Max minor axis
    h_mass_pbp = {}; % hydrometeor calculated mass
    h_density = {}; % hydrometeor density
    all_h_appears_ind = []; % initial time index of each snowflake in time array

    % Loop over all Hydrometeors for each frame: 
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
        % integration ~ need to double check with google search: 
        propstemp = regionprops(h_area_tmp_bool, 'PixelIdxList');
        % loop over some region property definded in previous line: 
            for jj = 1:numel(propstemp)
                h_max_area{end+1} = max(h_area_tmp(propstemp(jj).PixelIdxList));
                h_max_circ_area{end+1} = max(h_ellipse_area_tmp(propstemp(jj).PixelIdxList));
                h_max_maj_axis{end+1} = max(h_maj_axis_tmp(propstemp(jj).PixelIdxList));
                h_max_min_axis{end+1} = max(h_min_axis_tmp(propstemp(jj).PixelIdxList));
                h_delta_time{end+1} = numel(propstemp(jj).PixelIdxList);
                h_delta_temp_max{end+1} = max(h_dT_tmp(propstemp(jj).PixelIdxList));
                h_delta_temp_mean{end+1} = mean(h_dT_tmp(propstemp(jj).PixelIdxList));
                % multiply the area by the temperature difference for that hydrometeor:
                h_area_times_dT_pbp = h_area_tmp(propstemp(jj).PixelIdxList) .* h_dT_tmp(propstemp(jj).PixelIdxList);
                % now integrate:
                h_mass_pbp{end+1} = sum(k_dLv * h_area_times_dT_pbp);
            end
        end
    end
    
    %% Convert to matrixes:
    h_mass_pbp=cell2mat(h_mass_pbp); % hydrometeor mass
    h_max_area=cell2mat(h_max_area); % actual area
    h_max_circ_area=cell2mat(h_max_circ_area); % ellipse area
    h_max_maj_axis=cell2mat(h_max_maj_axis); % height
    h_max_min_axis=cell2mat(h_max_min_axis); % width
    h_delta_temp_max=cell2mat(h_delta_temp_max); % temperature diffrence between plate and water droplet using max intensity
    h_delta_temp_mean=cell2mat(h_delta_temp_mean); % temperature diffrence between plate and water droplet using mean intensity 
    % now multiples mass by "dt":
    h_mass_pbp = h_mass_pbp / vid_fps; 
    h_diam = (4/pi) * h_max_area.^(1/2); % convert area to diameter of hydrometeor
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

    %% organizes data into tables used for both output and post processing
    
    % handles FBF SWE data (used to compute SWE-factor later on)
    % fbf_table_raw = table(time_series_fbf', SWE_fbf);
    % fbf_table_raw.Properties.VariableNames{'Var1'} = 'Time';
    % fbf_table_raw.Time = datetime(fbf_table_raw.Time);
    % fbf_table_raw = table2timetable(fbf_table_raw);

    % conversions and calculations using PBP data:
    terminal_vel = terminalVelocity(h_max_area, h_max_circ_area, h_mass_pbp); % terminal velocity 
    void_space = (h_max_area - h_max_circ_area).^(1/2); % void space      
    for i = 1:length(void_space)
        if void_space(i) > 0
            void_space(i) = void_space(i);
        else
            void_space(i) = 0;
        end
    end
    surface_area_eq = ((9*pi)/16)^(1/3) * (h_mass_pbp ./ rho_water).^(2/3); % surface area of equivalent water droplet (See POF Singh et al. 2023)
    complexity = h_max_circ_area ./ h_max_area; % complexity (See CRST Morrison et al. 2023)
    sdi = h_max_area ./ surface_area_eq ; % SDI (See CRST Morrison et al. 2023)
    SWE_pbp = c1 * h_mass_pbp ./ (rho_water * hp_area(file_i)); % [mm]
    SWE_factor = sum(SWE_fbf) / sum(SWE_pbp); % swe factor is used to adjust pbp SWE 
    SWE_fbf_particles = SWE_pbp * SWE_factor; % adjusted SWE_pbp is effectivley SWE_fbf
    SWE_factor_particles = SWE_pbp ./ SWE_fbf_particles; % now calculate SWE factor for all particles
    SWE_pbp_accumulated = cumsum(SWE_pbp); % accumulated SWE_pbp
    SWE_fbf_accumulated = cumsum(SWE_fbf_particles); % accumulated SWE_fbf    
    snow_pbp = rho_water * (SWE_pbp ./ h_rho_hfd); % [mm]
    snow_pbp_accumulated = cumsum(snow_pbp); % [mm]
    snow_fbf = rho_water * (SWE_fbf_particles ./ h_rho_hfd); % [mm]
    snow_fbf_accumulated = cumsum(snow_fbf); % [mm]

    % store PBP data as table 
    pbp_table_particles = table(h_initial_time', h_evap_time', ... 
        h_mass_pbp', h_diam',  h_max_area', h_max_circ_area', ... 
        h_rho_sph', h_rho_hfd', h_vol_hfd', h_vol_sph', ... 
        h_delta_temp_max', h_delta_temp_mean', complexity', ... 
        sdi', SWE_pbp', SWE_fbf_particles', SWE_pbp_accumulated', ... 
        SWE_fbf_accumulated', snow_pbp', snow_fbf', snow_pbp_accumulated', ... 
        snow_fbf_accumulated', SWE_factor_particles');
    pbp_table_particles.Properties.VariableNames = {'Time', 'Evap Time (s)', 'Mass (kg)', 'Eff Diameter (m)',  ...
        'Max Area (m^2)', 'Max Circ Area (m^2)', 'Spherical Density (kg/m^3)', 'Heat Flux Density (kg/m^3)', ... 
        'Heat Flux Volume (m^3)', 'Spherical Volume (m^3)', 'Delta Temp Max', 'Delta Temp mean', ...
        'Complexity','SDI','PBP SWE (mm)','FBF SWE(mm)', 'PBP SWE Accumlation (mm)', ...
        'FBF SWE Accumulation (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'PBP Snow Accumulation (mm), ...' ...
        'FBF Snow Accumulation (mm)', 'SWE factor'};
    pbp_table_particles = table2timetable(pbp_table_particles); % convert to timetable 
    pbp_table_particles = sortrows(pbp_table_particles, 'Time'); % sort by time 
    pbp_table_particles_cell{file_i} = pbp_table_particles; % store as a cell so can be used outside of parforloop

 end

% stored everything as a cell to access out of parfor loop; now put back
% into one large time table 

% filters data to find where 0 < mass < .005 and 1/15 s < evaporation time < 60 s to omit residue on plate
% g1 = find((pbp_table_particles_cell.mass > 0 & ...
%     pbp_table_particles_cell.mass < residue_filter) & ...
%     (pbp_table_particles_cell.evap_time > evapTime_min & ...
%     pbp_table_particles_cell.evap_time < evapTime_max));

% if swe factor is abnormally high, adjust it to account for residue:
% if swe_factor > 1.95
%     swe_factor = 1.95; 
% end

% start_end_time_table = table(file_names', vid_start_time, vid_end_time); 
% start_end_time_table.length = start_end_time_table.vid_end_time - start_end_time_table.vid_start_time; 
% start_end_time_table.Properties.VariableNames{1} = 'file_name';