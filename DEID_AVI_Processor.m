%% DEID AVI File Processing Code
% Outputs particle by particle csv file to be post processed by
% DEID_AVI_Processor.m
% AUTHOR : Dhiraj Singh, Benjamin Silberman, Travis Morrison, Alex Blackmer
% LAST UPDATED: 04/10/2024
clear, clc, close all

%% Set path, obtain .avi file:
workingDir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID/';
inputDir = 'Atwater/JAN/Jan_09_10_storm/';
outputDir = 'processed_output_data/';
% Sets path for .m to run in:
cd([workingDir, inputDir]) 
% Check if the output directory exists
if ~exist(outputDir, 'dir')
    % If the directory doesn't exist, create it
    mkdir(outputDir);
    disp(['Directory "', outputDir, '" created successfully.']);
else
    % If the directory already exists, display a message
    disp(['Directory "', outputDir, '" already exists.']);
end

%% Identifies all .avi files in inputDir 
% Get a list of all files and folders in this folder
filesAndFolders = dir(".");
% Initialize an empty cell array to store the file names
fileNames = {};
% Loop through each item in the input directory
for file_i = 1:length(filesAndFolders)
    % Check if the item is a file (not a folder) and if it ends with .avi
    if ~filesAndFolders(file_i).isdir && endsWith(filesAndFolders(file_i).name, '.avi', 'IgnoreCase', true)
        % Get the name of the file and append it to the list
        fileNames{end+1} = filesAndFolders(file_i).name;
    end
end

%% Set global varabiles and constants:
pix_to_m_conversion = 1.0080625e-7; % Some factor to map pixel to real space
Colorbar_image_indexes = [1 1 384 288]; % Location of colorbar in pixel locations
Colorbar_Kapton_image_indexes = [1 27 384-1 288-27]; % Location of Kapton tape in pixel locations
Minmum_thres = 70; % Minimum threshold number in image accepted rbg ([0 255]) scale
Min_Hydrometeor_Size = 10; % Minumum Hydrometeor size in pixels 
Colorbar_Max_Temp = 145; % Max temperature set in colorbar on the physical screen of the tir software
Sort_Threshold = 20; % This it the RMS threshold between succesive images of snowflakes used in the sortPostitions_v2- I really don't understand this
Minimum_Drop_Life = 0; % Minimum number of frames a drop has to be visable to be processed
k_dLv = 0.0029; % Calibration constant, in paper thermal conductivity (k) of water See sect 4.1 in Dihiraj's paper -> (k/d(_eff))/Latent heat of vaporazation [units?]
camera_fps = 15;% This is the frames per second (fps) of the camera
L_constant = 2.594*10^6; % Latent heat of vaporazation of water, should be a function of tempertaure (Look at Stull textbook) [J/kg]
Rho_Water = 1000; % kg / m^3
HeatFlux_Density_coeff = 6.4418e04; % constant obtained from lab denisty of ice 

%% Begin DEID processing:
for file_i = 1:length(fileNames)
    disp(['Processing File: ', fileNames{file_i}])
    vid=VideoReader(fileNames{file_i});
    
    %% Creates datetime timeseries for output file
    % Get end time of avi file
    video_dir = dir(fileNames{file_i});
    video_end_time = datetime([video_dir.date]);
    % Get the duration of the video in seconds
    videoLengthSeconds = vid.Duration;
    % Specify the frequency (1/15 of a second)
    frequency = seconds(1/camera_fps);
    % Calculate the number of time steps
    num_time_steps = floor(videoLengthSeconds * camera_fps); % 1/15 of a second
    % Create a time series of date times starting from (video end time - video duration) and ending at the video end time
    time_series = datetime(video_end_time - (0:num_time_steps) * frequency, 'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
    time_series = flip(time_series);  % Flips time series so it ends at the end time
    % Get number of frames, preallocate cell for data, and obtain date
    % information:
    Num_frames = vid.NumFrames;                               
    Hydrometeor_Data = cell(Num_frames,1); 
    
    %% "Frame by frame method"; this is how Dhiraj is obtaining SWE for each
    % .avi file 
    % Preallocate variables saved in loop for speed:
    plate_temperature = nan(Num_frames,1);
    Sum_Hydrometeor_Area_times_DeltaT = nan(Num_frames,1);
    % Enter loop to process images: 
    for frame_ii = 1:Num_frames
        % Get image:    
        current_frame = read(vid, frame_ii);  
        % Clean image to get plate temperature: 
        current_frame_gray = rgb2gray(current_frame); % Convert frame of interest to gray scale
        current_frame_gray_cropped_wKapton = imcrop(current_frame_gray,Colorbar_image_indexes);% Crop out colorbar
        plate_temperature(frame_ii) = max(max(double(current_frame_gray_cropped_wKapton))); % Dhiraj assumes max temperature in image is the plate temperature with Kapton tape
        % Clean orginal image: 
        current_frame_gray_cropped = imcrop(current_frame_gray, Colorbar_Kapton_image_indexes); % Back to orginal grayscale image... now remove colorbar and kapton tape from image
        current_frame_filtered = current_frame_gray_cropped > 70; % Removed below 70, on rbg ([0, 255]) scale 
        current_frame_filtered_filled = imfill(current_frame_filtered, 'Holes'); % Clean up Hydrometeors
        % Get Hydrometeor properties: 
        Hydrometeor_geo_properties = regionprops(current_frame_filtered_filled, 'Centroid', 'Area','BoundingBox'); % Returns the centroid, the area of each blob, and the bounding box (left, top, width, height).
        % If no properties are found, go to next frame: 
        if (isempty(Hydrometeor_geo_properties))
            continue;
        end
        % Build Hydrometeor propertie matrices from regionprops values: 
        Hydrometeor_Bounding_box = cat(1,Hydrometeor_geo_properties.BoundingBox); % Concat all values to build matrix of Bounding Box indices
        Hydrometeor_Centroid = round(cat(1, Hydrometeor_geo_properties.Centroid)); % Concat all values to build matrix of Centroids
        Hydrometeor_Area_pix = cat(1, Hydrometeor_geo_properties.Area); % Concat all values to build matrix of Hydrometeor areas in pixels
        % Convert Hydrometeor areas to m^2: 
        Hydrometeor_Area = Hydrometeor_Area_pix.*pix_to_m_conversion; 
        Hydrometeor_Major_Axis = Hydrometeor_Bounding_box(:, 3) * sqrt(pix_to_m_conversion); % Hydrometeor major axis
        Hydrometeor_Minor_Axis = Hydrometeor_Bounding_box(:, 4) * sqrt(pix_to_m_conversion); % Hydrometeor minor axis
        % Get the Hydrometeor ellipse area:
        Hydrometeor_ellipse_area = Hydrometeor_Major_Axis.*Hydrometeor_Minor_Axis;
        % Calculating the difference in temperature of each centroid and
        % the plate:
        Hydrometeor_Centroid_index = sub2ind(size(current_frame_gray_cropped), Hydrometeor_Centroid(:, 2), Hydrometeor_Centroid(:, 1)); % Find the index of the centriods in orginal image
        Hydrometeor_Centroid_values = double(current_frame_gray_cropped(Hydrometeor_Centroid_index)); % Intensities of Centroid pixels of snow
        Plate_Hydrometeor_DeltaT = Colorbar_Max_Temp-(Hydrometeor_Centroid_values.*(Colorbar_Max_Temp/plate_temperature(frame_ii)));
        % Now calculate product of Hydrometeor area with the temp
        % difference:
        Hydrometeor_Area_times_DeltaT = Hydrometeor_Area.*Plate_Hydrometeor_DeltaT;         
        % Sum the product of individual area and temp. diff in each frame:
        % THIS IS HOW WE GET Hydrometeor MASS (FRAME BY FRAME METHOD)
        Sum_Hydrometeor_Area_times_DeltaT(frame_ii)=sum(Hydrometeor_Area_times_DeltaT);         
        % Build large matrix of Hydrometeor data
        Hydrometeor_Data{frame_ii} = cat(2, Hydrometeor_Centroid, Hydrometeor_Area, Plate_Hydrometeor_DeltaT, Hydrometeor_ellipse_area, Hydrometeor_Major_Axis, Hydrometeor_Minor_Axis); 
    end
    
    %% Frame by frame SWE calculation
    % Use current_frame_gray_cropped to calculate the area of hot plate:
    HotPlate_Area = size(current_frame_gray_cropped,1) * size(current_frame_gray_cropped,2) * pix_to_m_conversion;
    Hydrometeor_Mass_fbf = (k_dLv*Sum_Hydrometeor_Area_times_DeltaT) / camera_fps; % mass evaporates in each frame
    SWE_fbf = Hydrometeor_Mass_fbf / HotPlate_Area;
    SWE_fbf_accumulation = sum(SWE_fbf); 
    time_series_fbf = time_series(1:length(SWE_fbf));
    
    %% Sorting one frame to others frame ~ Data cleaning of some sorts: 
    % Create new cell array for sorted data: 
    Hydrometeor_Data_Sorted = Hydrometeor_Data; 
    % Loop over every data cell corresponding to each frame:
    % Uses "sortPositions_v2.m" - a code that Dhiraj wrote 
    for frame_jj = 2:Num_frames
        Hydrometeor_Data_Sorted{frame_jj} = sortPositions_v2(Hydrometeor_Data_Sorted{frame_jj-1}, Hydrometeor_Data_Sorted{frame_jj}, Sort_Threshold);
    end
    % Return frame with max number of Hydrometeors 
    Max_Hydrometeors_Obs = max(cellfun(@(x) size(x, 1), Hydrometeor_Data_Sorted, 'UniformOutput', 1));
    % Now pad the data with zeros so the frames all have the same number:  
    Hydrometeor_Data_Sorted = cellfun(@(x) cat(1, x, zeros(Max_Hydrometeors_Obs-size(x, 1), 7)),...
        Hydrometeor_Data_Sorted, 'UniformOutput', 0);
    
    %%  Isolating the variables and put them into a matrix to work with:
    % For reference: [Hydrometeor_Centroid, Hydrometeor_Area,
    % Plate_Hydrometeor_DeltaT, Hydrometeor_ellipse_area,
    % Hydrometeor_Major_Axis,Hydrometeor_Minor_Axis]
    Final_Hydrometeor_Area = cellfun(@(x) x(:, 3), Hydrometeor_Data_Sorted, 'UniformOutput', 0);
    Final_Delta_T = cellfun(@(x) x(:, 4), Hydrometeor_Data_Sorted, 'UniformOutput', 0);
    Final_ellipse_Area = cellfun(@(x) x(:, 5), Hydrometeor_Data_Sorted, 'UniformOutput', 0);
    Final_Major_Axis = cellfun(@(x) x(:, 6), Hydrometeor_Data_Sorted, 'UniformOutput', 0);
    Final_Minor_Axis = cellfun(@(x) x(:, 7), Hydrometeor_Data_Sorted, 'UniformOutput', 0);
    %Convert to matrix: 
    Final_Hydrometeor_Area = cat(2,Final_Hydrometeor_Area{:});
    Final_Delta_T = cat(2,Final_Delta_T{:});
    Final_ellipse_Area = cat(2,Final_ellipse_Area{:});
    Final_Major_Axis = cat(2,Final_Major_Axis{:});
    Final_Minor_Axis = cat(2,Final_Minor_Axis{:});

    %% "Particle by particle method"
    % Calulcate the evolution of each Hydrometeor to get properties.
    % This method 'cleans' the plate to not double count Hydrometeors that
    % appear across multiple frames. 
    % Preallocate variables:
    Hydrometeor_Delta_time = {}; % Hydrometeor time to melt and evaporate
    Hydrometeor_Delta_T = {}; % Hydrometeor centroid's max intensity
    Hydrometeor_Delta_T1 = {}; % Hydrometeor centroid's mean intensity
    Hydrometeor_Max_Area = {}; % Hydrometeor maximum area
    Hydrometeor_Max_Circumscribed_Area = {};% Hydrometeor Max circumscribed area
    Hydrometeor_Max_Major_Axis = {};% Hydrometeor Max Major axis
    Hydrometeor_Max_Minor_Axis = {};% Hydrometeor Max minor axis
    Hydrometeor_Mass_pbp = {}; % Hydrometeor calculated mass
    Hydrometeor_Density = {}; % Hydrometeor density
    Hydrometeor_inital_time_index=[]; % Index of each snowflake in time array

    %% Loop over all Hydrometeors: 
    % This is where the cleaning occurs!
    for Hydrometeor_ii = 1:Max_Hydrometeors_Obs
        Hydrometeor_appear_evap_boolean = diff(Final_Hydrometeor_Area(Hydrometeor_ii, :) > 0); 
        % Create array of indexes which are booleans, +1 is when snowflake starts, and -1 when it leaves
        % diff is [X(3)-X(2)] for backwward FD
        if sum(Hydrometeor_appear_evap_boolean)<0
            continue;
        end
        Hydrometeor_appears_index = find(Hydrometeor_appear_evap_boolean > 0); % finds the index when the Hydrometeor appears
        Hydrometeor_evaps_index = find(Hydrometeor_appear_evap_boolean < 0); % Hydrometeor dissapears index
        Hydrometeor_inital_time_index=cat(2,Hydrometeor_inital_time_index,Hydrometeor_appears_index); % time when snow arrives the plate
        if isempty(Hydrometeor_evaps_index) % get rid of Hydrometeors which do no evaporate completely
            continue; 
        else
        % Isolate the Hydrometeor properties for the lifetime of the
        % Hydrometeor: 
        Hydrometeor_Delta_T_tmp = Final_Delta_T(Hydrometeor_ii, Hydrometeor_appears_index(1):Hydrometeor_evaps_index(end)+1);
        Hydrometeor_ellipse_Area_tmp = Final_ellipse_Area(Hydrometeor_ii, Hydrometeor_appears_index(1):Hydrometeor_evaps_index(end)+1);
        Hydrometeor_Major_Axis_tmp = Final_Major_Axis(Hydrometeor_ii,Hydrometeor_appears_index(1):Hydrometeor_evaps_index(end)+1);
        Hydrometeor_Minor_Axis_tmp = Final_Minor_Axis(Hydrometeor_ii, Hydrometeor_appears_index(1):Hydrometeor_evaps_index(end)+1);
        % Prepare area of each Hydrometeor for temporal intergration: 
        Hydrometeor_area_tmp = Final_Hydrometeor_Area(Hydrometeor_ii, Hydrometeor_appears_index(1):Hydrometeor_evaps_index(end)+1); 
        Hydrometeor_area_tmp_bool = Hydrometeor_area_tmp > 0; % Isolate for positive values
        Hydrometeor_area_tmp_bool = bwareaopen(Hydrometeor_area_tmp_bool, Minimum_Drop_Life); % Any Hydrometeor which lives for less than 'Minimum_Drop_Life' is discarded
        % Integration ~ need to double check with google search: 
        propstemp = regionprops(Hydrometeor_area_tmp_bool, 'PixelIdxList');
        % Loop over some region property definded in previous line: 
            for jj = 1:numel(propstemp)
                Hydrometeor_Max_Area{end+1} = max(Hydrometeor_area_tmp(propstemp(jj).PixelIdxList));
                Hydrometeor_Max_Circumscribed_Area{end+1} = max(Hydrometeor_ellipse_Area_tmp(propstemp(jj).PixelIdxList));
                Hydrometeor_Max_Major_Axis{end+1} = max(Hydrometeor_Major_Axis_tmp(propstemp(jj).PixelIdxList));
                Hydrometeor_Max_Minor_Axis{end+1} = max(Hydrometeor_Minor_Axis_tmp(propstemp(jj).PixelIdxList));
                Hydrometeor_Delta_time{end+1} = numel(propstemp(jj).PixelIdxList);
                Hydrometeor_Delta_T{end+1} = max(Hydrometeor_Delta_T_tmp(propstemp(jj).PixelIdxList));
                Hydrometeor_Delta_T1{end+1} = mean(Hydrometeor_Delta_T_tmp(propstemp(jj).PixelIdxList));
                % Multiply the area by the temperature difference for that Hydrometeor:
                Hydrometeor_Area_times_DeltaT_pbp = Hydrometeor_area_tmp(propstemp(jj).PixelIdxList).*Hydrometeor_Delta_T_tmp(propstemp(jj).PixelIdxList);
                % Now integrate:
                Hydrometeor_Mass_pbp{end+1} = trapz(k_dLv*Hydrometeor_Area_times_DeltaT_pbp);
            end
        end
    end
    
    %% Convert to matrixes:
    Hydrometeor_Mass_pbp=cell2mat(Hydrometeor_Mass_pbp); % Hydrometeor mass
    Hydrometeor_Max_Area=cell2mat(Hydrometeor_Max_Area); % Actual area
    Hydrometeor_Max_Circumscribed_Area=cell2mat(Hydrometeor_Max_Circumscribed_Area); % Ellipse area
    Hydrometeor_Max_Major_Axis=cell2mat(Hydrometeor_Max_Major_Axis); % Height
    Hydrometeor_Max_Minor_Axis=cell2mat(Hydrometeor_Max_Minor_Axis); % Width
    Hydrometeor_Delta_T=cell2mat(Hydrometeor_Delta_T); % Temperature diffrence between plate and water droplet using max intensity
    Hydrometeor_Delta_T1=cell2mat(Hydrometeor_Delta_T1); % Temperature diffrence between plate and water droplet using mean intensity 
    % Now multiples mass by "dt":
    Hydrometeor_Mass_pbp = Hydrometeor_Mass_pbp/camera_fps; 
    Hydrometeor_Diameter = (1.12)*Hydrometeor_Max_Area.^0.5; % Convert area to diameter of Hydrometeor
    Water_eq_Diameter = (6.*Hydrometeor_Mass_pbp/3140).^0.33; % Water equi. diameter. Not sure where this is from
    Time_2_evap = cell2mat(Hydrometeor_Delta_time)*(1/camera_fps); % Evaporation time
    Time_2_evap_norm = sort(Hydrometeor_inital_time_index)/camera_fps; % Time in sec to what? Evap?
    Spherical_Volume = 0.75*Hydrometeor_Max_Area.^1.5; % Spherical volume
    Sphere_Density = Hydrometeor_Mass_pbp./Spherical_Volume; % Density calculation: spherical assumption
    Energy_perTime = Hydrometeor_Mass_pbp*L_constant./(Hydrometeor_Max_Area.*Time_2_evap); % Heat flux method: energy per unit area per time
    HeatFlux_Density1 = HeatFlux_Density_coeff * Hydrometeor_Mass_pbp./(Hydrometeor_Max_Area.*Time_2_evap.*Hydrometeor_Delta_T);
    HeatFlux_Density2 = HeatFlux_Density_coeff * Hydrometeor_Mass_pbp./(Hydrometeor_Max_Area.*Time_2_evap.*Hydrometeor_Delta_T1);
    HeatFlux_Volume = Hydrometeor_Mass_pbp./HeatFlux_Density2; % Volume of each snowflakes using heat flux method density
    Hydrometeor_initial_time_indexes= Hydrometeor_inital_time_index(1:length(Hydrometeor_Mass_pbp));
    Hydrometeor_initial_times = time_series(Hydrometeor_initial_time_indexes);
    Hydrometeor_Time_2_evap_norm = Time_2_evap_norm(1:length(Hydrometeor_Mass_pbp));

    %% Particle by particle SWE calculation
    SWE_pbp = Hydrometeor_Mass_pbp ./ HotPlate_Area;
    SWE_pbp_accumulation = sum(SWE_pbp);
    SWE_factor = SWE_fbf_accumulation / SWE_pbp_accumulation; % we use this factor when using SWE to calculate Snow Depth
    mean_density = mean(HeatFlux_Density1, 'omitnan'); % Desnity-Heat Flux Method
    mean_density1 = mean(HeatFlux_Density2, 'omitnan'); % Desnity-Heat Flux Method
    
    %% Organizes data into output arrays
    Data_FBF= table(time_series_fbf', SWE_fbf);
    Data_PBP = table(Hydrometeor_initial_times', Hydrometeor_initial_time_indexes',Hydrometeor_Time_2_evap_norm', Hydrometeor_Mass_pbp', Hydrometeor_Diameter',  ...
        Hydrometeor_Max_Area', Hydrometeor_Max_Circumscribed_Area', Time_2_evap', Sphere_Density', HeatFlux_Density1', HeatFlux_Density2', Energy_perTime', ... 
        HeatFlux_Volume', Water_eq_Diameter', Hydrometeor_Max_Major_Axis', Hydrometeor_Max_Minor_Axis',  Hydrometeor_Delta_T', Hydrometeor_Delta_T1');
    Data_PBP.Properties.VariableNames = {'Hydrometeor_initial_times', 'Hydrometeor_initial_time_indexes','Hydrometeor_Time_2_evap_norm', 'Hydrometeor_Mass_pbp', 'Hydrometeor_Diameter',  ...
        'Hydrometeor_Max_Area', 'Hydrometeor_Max_Circumscribed_Area', 'Time_2_evap', 'Sphere_Density', 'HeatFlux_Density1', 'HeatFlux_Density2', 'Energy_perTime', ... 
        'HeatFlux_Volume', 'Water_eq_Diameter', 'Hydrometeor_Max_Major_Axis', 'Hydrometeor_Max_Minor_Axis',  'Hydrometeor_Delta_T', 'Hydrometeor_Delta_T1'};
    % Data_PBP.Properties.VariableNames = {'Initial_Time', 'Initial_Time_Index','Time_2_Evap_Norm', 'Mass', 'Diameter',  ...
    %     'Max_Area', 'Max_Circumscribed_Area', 'Time_2_Evap', 'Sphere_Density', 'Heat_Flux_Density_1', 'Heat_Flux_Density_2', 'Energy_Per_Time', ... 
    %     'Heat_Flux_Volume', 'Water_Eq_Diameter', 'Max_Major_Axis', 'Max_Minor_Axis',  'Delta_T', 'Delta_T1'};

    %% Saves processed output arrays to .csv for single .avi file
    % Replace the file extension with '.csv'
    outputFileName = strrep(fileNames{file_i}, '.avi', '');
    writetable(Data_PBP, [outputDir, outputFileName, '_PBP.csv']);
    writetable(Data_FBF, [outputDir, outputFileName, '_FBF.csv']);
    disp(['Saved Output for: ', outputFileName])
end


