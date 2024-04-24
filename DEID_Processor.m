%% DEID AVI File Processing Code
% Outputs particle by particle csv file to be post processed by
% DEID_AVI_Processor.m
% AUTHOR : Dhiraj Singh, Benjamin Silberman, Travis Morrison, Alex Blackmer
% LAST UPDATED: 04/10/2024
%                                                                     
clear, clc, close all


%% Set Up 
% Preallocate Storage
col_names = {'timestamp', 'hfd_1', 'hfd_2', 'SWE_FBF_mm', 'SWE_PBP_mm', 'SWE_PBP_F_mm', 'rho_hfd1_F', 'snow_PBP_hfd1_mm'};
col_types = {'datetime', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};

output_table = table('Size', [0, length(col_names)], ...
                         'VariableNames', col_names, ...
                         'VariableTypes', col_types);

% Set path, obtain .avi file:
working_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID/Atwater/JAN/test';
% working_dir = 'Z:\DEID\Atwater\JAN\test';
% Sets path for .m to run in:
cd(working_dir) 

%% Identifies all .avi files in inputDir 
% Get a list of all files and folders in this folder
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

%% Set global varabiles and constants:
% Specifies resampling period frequency
time_step = minutes(5);      
% Defines global variables
mm_to_inches = 1/25.4;
pix_to_m_conversion = 1.0080625e-7; % Some factor to map pixel to real space
colorbar_image_indexes = [1 1 384 288]; % Location of colorbar in pixel locations
colorbar_kapton_image_indexes = [1 27 384-1 288-27]; % Location of Kapton tape in pixel locations
colorbar_max_temp = 145; % Max temperature set in colorbar on the physical screen of the tir software
min_thres = 70; % Minimum threshold number in image accepted rbg ([0 255]) scale
sort_threshold = 20; % This it the RMS threshold between succesive images of snowflakes used in the sortPostitions_v2- I really don't understand this
min_h_size = 10; % Minumum Hydrometeor size in pixels 
minimum_drop_life = 0; % Minimum number of frames a drop has to be visable to be processed
k_dLv = 0.0029; % Calibration constant, in paper thermal conductivity (k) of water See sect 4.1 in Dihiraj's paper -> (k/d(_eff))/Latent heat of vaporazation [units?]
l_constant = 2.594*10^6; % Latent heat of vaporazation of water, should be a function of tempertaure (Look at Stull textbook) [J/kg]
rho_water = 1000; % Density of water [kg/m^3]
hf_rho_coeff = 6.4418e04; % Heat flux density constant obtained from lab denisty of ice. [(K*s)/m]

%% Begin DEID processing:
for file_i = 1:length(file_names)
    disp(['Processing File: ', file_names{file_i}])
    vid=VideoReader(file_names{file_i});
    
    %% Creates datetime timeseries for output file
    % Get end time of avi file
    vid_dir = dir(file_names{file_i});
    vid_end_time = datetime([vid_dir.date]);
    % Get the duration of the video in seconds
    vid_length = vid.Duration;
    vid_fps = vid.FrameRate;
    % Specify the frequency 
    frequency = seconds(1/vid_fps);
    % Calculate the number of time steps
    num_frames = vid.NumFrames;
    % Create a time series of date times starting from (video end time - video duration) and ending at the video end time
    time_series = datetime(vid_end_time - (0:num_frames) * frequency, 'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
    time_series = flip(time_series);  % Flips time series so it ends at the end time
    % Get number of frames, preallocate cell for data, and obtain date information:
    h_data = cell(num_frames,1);    % Hydrometeor data
    
    %% "Frame by frame method"; this is how Dhiraj is obtaining SWE for each .avi file 
    % Preallocate variables saved in loop for speed:
    plate_temp = nan(num_frames,1);
    sum_h_area_times_dt = nan(num_frames,1);
    % Enter loop to process images: 
    for frame_ii = 1:num_frames 
        % Get image:    
        current_frame = read(vid, frame_ii);  
        % Clean image to get plate temperature: 
        current_frame_gray = rgb2gray(current_frame); % Convert frame of interest to gray scale
        current_frame_gray_cropped_wKapton = imcrop(current_frame_gray,colorbar_image_indexes);% Crop out colorbar
        plate_temp(frame_ii) = max(max(double(current_frame_gray_cropped_wKapton))); % Dhiraj assumes max temperature in image is the plate temperature with Kapton tape
        % Clean orginal image: 
        current_frame_gray_cropped = imcrop(current_frame_gray, colorbar_kapton_image_indexes); % Back to orginal grayscale image... now remove colorbar and kapton tape from image
        current_frame_filtered = current_frame_gray_cropped > min_thres; % Removed below min threshold, on rbg ([0, 255]) scale 
        current_frame_filtered_filled = imfill(current_frame_filtered, 'Holes'); % Clean up Hydrometeors
        % Get Hydrometeor properties: 
        h_geo_prop = regionprops(current_frame_filtered_filled, 'Centroid', 'Area','BoundingBox'); % Returns the centroid, the area of each blob, and the bounding box (left, top, width, height).
        % If no properties are found, go to next frame: 
        if (isempty(h_geo_prop))
            continue;
        end
        % Build Hydrometeor propertie matrices from regionprops values: 
        h_bounding_box = cat(1,h_geo_prop.BoundingBox); % Concat all values to build matrix of Bounding Box indices
        h_centroid = round(cat(1, h_geo_prop.Centroid)); % Concat all values to build matrix of Centroids
        h_area_pix = cat(1, h_geo_prop.Area); % Concat all values to build matrix of Hydrometeor areas in pixels
        % Convert Hydrometeor areas to m^2: 
        h_area = h_area_pix.*pix_to_m_conversion; 
        h_major_axis = h_bounding_box(:, 3) * sqrt(pix_to_m_conversion); % Hydrometeor major axis
        h_minor_axis = h_bounding_box(:, 4) * sqrt(pix_to_m_conversion); % Hydrometeor minor axis
        % Get the Hydrometeor ellipse area:
        h_elipse_area = h_major_axis.*h_minor_axis;
        % Calculating the difference in temperature of each centroid and
        % the plate:
        h_centroid_i = sub2ind(size(current_frame_gray_cropped), h_centroid(:, 2), h_centroid(:, 1)); % Find the index of the centriods in orginal image
        h_centroid_values = double(current_frame_gray_cropped(h_centroid_i)); % Intensities of Centroid pixels of snow
        plate_h_dtemp = colorbar_max_temp-(h_centroid_values.*(colorbar_max_temp/plate_temp(frame_ii)));
        % Now calculate product of Hydrometeor area with the temp
        % difference:
        h_area_times_dtemp = h_area.*plate_h_dtemp;         
        % Sum the product of individual area and temp. diff in each frame:
        % THIS IS HOW WE GET Hydrometeor MASS (FRAME BY FRAME METHOD)
        sum_h_area_times_dt(frame_ii)=sum(h_area_times_dtemp);         
        % Build large matrix of Hydrometeor data
        h_data{frame_ii} = cat(2, h_centroid, h_area, plate_h_dtemp, h_elipse_area, h_major_axis, h_minor_axis); 
    end

    
    %% Frame by frame SWE calculation
    % Use current_frame_gray_cropped to calculate the area of hot plate:
    hp_area = size(current_frame_gray_cropped,1) * size(current_frame_gray_cropped,2) * pix_to_m_conversion;    % Hotplate Area [m^2]
    h_mass_fbf = (k_dLv*sum_h_area_times_dt) / vid_fps; % mass evaporates in each frame
    SWE_FBF_mm = h_mass_fbf / hp_area;
    SWE_fbf_accumulation = sum(SWE_FBF_mm); 
    time_series_fbf = time_series(1:length(SWE_FBF_mm));
    
    %% Sorting one frame to others frame ~ Data cleaning of some sorts: 
    % Create new cell array for sorted data: 
    h_data_sorted = h_data; 
    % Loop over every data cell corresponding to each frame:
    % Uses "sortPositions_v2.m" - a code that Dhiraj wrote 
    for frame_jj = 2:num_frames
        h_data_sorted{frame_jj} = sortPositions_v2(h_data_sorted{frame_jj-1}, h_data_sorted{frame_jj}, sort_threshold);
    end
    % Return frame with max number of Hydrometeors 
    max_h_obs = max(cellfun(@(x) size(x, 1), h_data_sorted, 'UniformOutput', 1));
    % Now pad the data with zeros so the frames all have the same number:  
    h_data_sorted = cellfun(@(x) cat(1, x, zeros(max_h_obs-size(x, 1), 7)),...
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

    %% Loop over all Hydrometeors: 
    % This is where the cleaning occurs!
    for h_ii = 1:max_h_obs
        h_appear_evap_bool = diff(h_area_final(h_ii, :) > 0); 
        % Create array of indexes which are booleans, +1 is when snowflake starts, and -1 when it leaves
        % diff is [X(3)-X(2)] for backwward FD
        if sum(h_appear_evap_bool)<0
            continue;
        end
        h_appears_ind = find(h_appear_evap_bool > 0); % finds the index when the Hydrometeor appears
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
                h_area_times_dT_pbp = h_area_tmp(propstemp(jj).PixelIdxList).*h_dT_tmp(propstemp(jj).PixelIdxList);
                % Now integrate:
                h_mass_pbp{end+1} = trapz(k_dLv*h_area_times_dT_pbp);
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
    h_diam = (1.12) * h_max_area.^0.5; % Convert area to diameter of Hydrometeor
    water_eq_diameter = (6 .* h_mass_pbp / 3140).^0.33; % Water equi. diameter. Not sure where this is from
    time_2_evap = cell2mat(h_delta_time) * (1/vid_fps); % Evaporation time
    time_2_evap_norm = sort(h_init_time_ind) / vid_fps; % Time in sec to what? Evap?
    spherical_volume = 0.75 * h_max_area.^1.5; % Spherical volume
    sphere_density = h_mass_pbp ./ spherical_volume; % Density calculation: spherical assumption
    energy_per_time = h_mass_pbp * l_constant ./ (h_max_area.*time_2_evap); % Heat flux method: energy per unit area per time
    heat_flux_density_1 = hf_rho_coeff * h_mass_pbp ./ (h_max_area .* time_2_evap .* h_delta_temp_max);
    heat_flux_density_2 = hf_rho_coeff * h_mass_pbp ./ (h_max_area .* time_2_evap .* h_delta_temp_mean);
    heat_flux_volume = h_mass_pbp ./ heat_flux_density_2; % Volume of each snowflakes using heat flux method density
    hydrometeor_initial_time_indexes = h_init_time_ind(1:length(h_mass_pbp));
    hydrometeor_initial_times = time_series(hydrometeor_initial_time_indexes);
    hydrometeor_Time_2_evap_norm = time_2_evap_norm(1:length(h_mass_pbp));

    %% Particle by particle SWE calculation
    SWE_pbp = h_mass_pbp ./ hp_area;
    SWE_pbp_accumulation = sum(SWE_pbp);
    SWE_factor = SWE_fbf_accumulation / SWE_pbp_accumulation; % we use this factor when using SWE to calculate Snow Depth
    
    %% Organizes data into output arrays
    Data_FBF= table(time_series_fbf', SWE_FBF_mm);
    Data_PBP = table(hydrometeor_initial_times', hydrometeor_initial_time_indexes',hydrometeor_Time_2_evap_norm', h_mass_pbp', h_diam',  ...
        h_max_area', h_max_circ_area', time_2_evap', sphere_density', heat_flux_density_1', heat_flux_density_2', energy_per_time', ... 
        heat_flux_volume', water_eq_diameter', h_max_maj_axis', h_max_min_axis',  h_delta_temp_max', h_delta_temp_mean');
    Data_PBP.Properties.VariableNames = {'Hydrometeor_initial_times', 'Hydrometeor_initial_time_indexes','Hydrometeor_Time_2_evap_norm', 'Hydrometeor_Mass_pbp', 'Hydrometeor_Diameter',  ...
        'Hydrometeor_Max_Area', 'Hydrometeor_Max_Circumscribed_Area', 'Time_2_evap', 'Sphere_Density', 'HeatFlux_Density1', 'HeatFlux_Density2', 'Energy_perTime', ... 
        'HeatFlux_Volume', 'Water_eq_Diameter', 'Hydrometeor_Max_Major_Axis', 'Hydrometeor_Max_Minor_Axis',  'Hydrometeor_Delta_T', 'Hydrometeor_Delta_T1'};
    % Data_PBP.Properties.VariableNames = {'Initial_Time', 'Initial_Time_Index','Time_2_Evap_Norm', 'Mass', 'Diameter',  ...
    %     'Max_Area', 'Max_Circumscribed_Area', 'Time_2_Evap', 'Sphere_Density', 'Heat_Flux_Density_1', 'Heat_Flux_Density_2', 'Energy_Per_Time', ... 
    %     'Heat_Flux_Volume', 'Water_Eq_Diameter', 'Max_Major_Axis', 'Max_Minor_Axis',  'Delta_T', 'Delta_T1'};



    %% Post Processing Script starts here
    % Handles FBF SWE data
    fbf_table_raw = Data_FBF;
    fbf_table_raw.Properties.VariableNames{'Var1'} = 'timestamp';
    fbf_table_raw.timestamp = datetime(fbf_table_raw.timestamp);
    fbf_table_raw = table2timetable(fbf_table_raw);
    % Resamples time series at desired interval
    fbf_table = retime(fbf_table_raw, 'regular', 'sum', 'TimeStep', time_step);
    fbf_table.SWE_FBF_accum_mm = cumsum(fbf_table.SWE_FBF_mm);
    
    %% Handles PBP Data
    % Reads in particle by particle table data
    pbp_table_raw = Data_PBP;
    % Sorts table by timestamp
    timestamp = datetime(pbp_table_raw.Hydrometeor_initial_times);
    pbp_table_raw = sortrows(pbp_table_raw, 'Hydrometeor_initial_times');
    % These variables can be replaced when this script is merged with
    % DEID_AVI_PROCESSOR.m
    mass = pbp_table_raw.Hydrometeor_Mass_pbp;    
    diameter = pbp_table_raw.Hydrometeor_Diameter;    
    max_area = pbp_table_raw.Hydrometeor_Max_Area;    
    max_circ_area = pbp_table_raw.Hydrometeor_Max_Circumscribed_Area; 
    time_to_evap = pbp_table_raw.Time_2_evap;        
    delta_T = pbp_table_raw.Hydrometeor_Delta_T;    
    delta_T1 = pbp_table_raw.Hydrometeor_Delta_T1; 
    density_sph = pbp_table_raw.Sphere_Density;    
    % Reads in heat flux density directly 
    hfd_1 = pbp_table_raw.HeatFlux_Density1;
    hfd_2 = pbp_table_raw.HeatFlux_Density2;
    
    % Computes terminal velocity
    % Constants
    rho_air = 0.9;
    g = 9.8;
    pi = 3.14;
    eta = 1.81*10^-5;
    % Derived values
    Area_ratio = max_circ_area./max_area;
    nu = eta/rho_air;
    x1 = 8*mass.*g*rho_air;
    x2 = pi*eta^2;
    x3 = Area_ratio.^0.25;
    X = x1.*x3./x2;
    p1 = (1+0.1519.*X.^0.5).^0.5;
    Re = 8.5*(p1-1).^2;
    v = Re.*eta.*(pi./max_circ_area).^0.5;
    v_t = v./(2.*rho_air);
    vvol = mass./rho_water;
    A_m = 1.2*vvol.^(2/3);
    cpx1 = max_circ_area./max_area;
    sdi = max_area./A_m;
    mu = 1.5*10^-5;
    v_st = (1/(18*mu))* hfd_1.*diameter.^2;
    
    % Stores derived data in table for easy access
    pbp_table_raw = table(timestamp, mass, density_sph, hfd_1, hfd_2, v_t, max_area, cpx1, sdi, diameter, max_circ_area, time_to_evap);
    pbp_table_raw = sortrows(pbp_table_raw, 'timestamp');
    % Filters data to find where 0 < mass < .005
    % Not sure if code is still necessary after rewriting time averaging method
    [g1,g2] = find(pbp_table_raw.mass > 0 & pbp_table_raw.mass < 0.005);
    pbp_table_raw = pbp_table_raw(g1,:);
    
    % Appends volume data to table 
    pbp_table_raw.VV1 = pbp_table_raw.mass./(pbp_table_raw.density_sph);    % Volume calculated using spherical density
    pbp_table_raw.VV2 = pbp_table_raw.mass./(pbp_table_raw.hfd_1);          % Volume calculated using heat flux density 1
    pbp_table_raw.VV3 = pbp_table_raw.mass./(pbp_table_raw.hfd_2);          % Volume calculated using heat flux density 2
    pbp_table_raw = table2timetable(pbp_table_raw);
    
    %% Computes averaged and summed PBP data
    % Averages PBP data 
    avg_cols = {'density_sph', 'hfd_1', 'hfd_2', 'v_t', 'cpx1', 'sdi', 'diameter'};
    avg_table = retime(pbp_table_raw(:, avg_cols), 'regular', 'mean', 'TimeStep', time_step);
    % Sums PBP data
    sum_cols = {'mass', 'max_area', 'max_circ_area', 'VV1', 'VV2', 'VV3'};
    sum_table = retime(pbp_table_raw(:, sum_cols), 'regular', 'sum', 'TimeStep', time_step);
    
    pbp_table = horzcat(avg_table, sum_table);
    
    %% Total SWE per averaging period PBP data
    pbp_table.SWE_PBP_mm = pbp_table.mass ./ (hp_area);      
    pbp_table.SWE_PBP_accum_mm = cumsum(pbp_table.SWE_PBP_mm);  
    % Finds difference factor between FBF SWE and PBP SWE and adjusts PBP SWE
    factor = fbf_table.SWE_FBF_accum_mm(end) / pbp_table.SWE_PBP_accum_mm(end);
    pbp_table.SWE_PBP_F_mm = pbp_table.SWE_PBP_mm * factor;                                   
    pbp_table.SWE_PBP_F_accum_mm = cumsum(pbp_table.SWE_PBP_F_mm);            
    
    %% Total Snow per averaging period PBP data
    % Adjusted density
    pbp_table.rho_spd = pbp_table.mass ./ pbp_table.VV1;    % Density from spherical density method [kg/m^3]
    pbp_table.rho_hfd1 = pbp_table.mass ./ pbp_table.VV2;   % Density from HFD1 density method [kg/m^3]
    pbp_table.rho_hfd2 = pbp_table.mass ./ pbp_table.VV3;   % Density from HFD2 density method [kg/m^3]

    % Finds difference factor between volume methods
    f_spd_hfd1 = mean(pbp_table.VV1) ./ mean(pbp_table.VV2);   % Adjustment factor between spd and HFD1 density methods
    pbp_table.VV2_F = f_spd_hfd1 * pbp_table.VV2;                % Adjusted volume for HFD1
    pbp_table.rho_hfd1_F = pbp_table.mass ./ pbp_table.VV2_F;  % Density from HFD1 adjusted for spd density method [kg/m^3]
    pbp_table.snow_PBP_hfd1_mm = rho_water * pbp_table.SWE_PBP_F_mm ./ pbp_table.rho_hfd1_F;   
    pbp_table.snow_PBP_hfd1_acc_mm = cumsum(pbp_table.snow_PBP_hfd1_mm);
   
    f_spd_hfd2 = mean(pbp_table.VV1) ./ mean(pbp_table.VV3);   % Adjustment factor between spd and HFD2 density methods
    pbp_table.VV3_F = f_spd_hfd2 * pbp_table.VV3;                % Adjusted volume for HFD2
    pbp_table.rho_hfd2_F = pbp_table.mass ./ pbp_table.VV3_F;  % Density from HFD2 adjusted for spd density method [kg/m^3]
    pbp_table.snow_PBP_hfd2_mm = rho_water * pbp_table.SWE_PBP_F_mm ./ pbp_table.rho_hfd2_F;   
    pbp_table.snow_PBP_hfd2_acc_mm = cumsum(pbp_table.snow_PBP_hfd2_mm);
    
    % Appends data for single video to output table
    output = synchronize(fbf_table, pbp_table);
    output = timetable2table(output);
    % Loop through each variable  and replace NaNs with zeros.
    % NaNs are showing up for some videos becaus the pbp and fbf tables
    % dont have the same number of time stamps. Replacing Nans is a quick
    % fix but should be further explored.
    output_names = output.Properties.VariableNames;
    for i = 2:length(output_names)
        output.(output_names{i})(isnan(output.(output_names{i}))) = 0;
    end

    output_table = [output_table; output(:, col_names)];
end

%% Sorts and cumulatively sums data for SWE and Snow totals
output_table = sortrows(output_table, 'timestamp');
output_table.SWE_FBF_acc_mm = cumsum(output_table.SWE_FBF_mm);
output_table.SWE_PBP_acc_mm = cumsum(output_table.SWE_PBP_mm);
output_table.SWE_PBP_F_acc_mm = cumsum(output_table.SWE_PBP_F_mm);
output_table.snow_PBP_hfd1_acc_mm = cumsum(output_table.snow_PBP_hfd1_mm);
output_table.snow_PBP_hfd1_acc_in = output_table.snow_PBP_hfd1_acc_mm * mm_to_inches;


%% Saves processed output arrays to .csv for single .avi file
% Replace the file extension with '.csv'
currentDir = pwd;
[~, parent_dir, ~] = fileparts(currentDir);
writetable(output_table, [parent_dir, '.csv']);
disp(['Saved Output for: ', parent_dir])


