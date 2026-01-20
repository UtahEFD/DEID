%% calibrates k/d coefficient using 20ul water droplets 

% droplets placed on the thermal hotplate at atwater snow study plot

% last calibration date: 20260107 

clear, clc
%% set filepath, global variables, and physical constants.

working_dir = '/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2025_2026/kd_testing';
output_dir = working_dir;

%% global variables and physical constants

% unit conversions:

pix_to_m_conversion = .01/40; % m per pix 
pix_to_m2_conversion = pix_to_m_conversion^2; % m^2 per pix^2

% thresholds & filters:

min_thres = 70; % minimum threshold number in image accepted rbg ([0 255]) scale
minimum_hydro_area = 100; % the minimum number of pixels a hydrometeor must contain to be analyzed 
colorbar_max_temp = 145; % max temperature set in colorbar on the physical screen of the tir software
sort_threshold = 20; % this is the RMS threshold between succesive images of snowflakes used in the sortPostitions_v2. dhiraj calibrated this in the lab. ß

% DEID specific parameters:

colorbar_image_indexes = [1 1 384 288]; % location of colorbar in pixel locations
crop_index = 43; % use this to specify indices to crop out kapton tape
colorbar_kapton_image_indexes = [1 (colorbar_image_indexes(2)+crop_index) 383 (colorbar_image_indexes(4)-crop_index)]; % location of Kapton tape in pixel locations

% to calibrate k/d, we know mass, so we instead solve the adjusted conservation of energy equation for k/d coefficient  
mass_calib = 2e-05; % 20uL of water 

%% move to working directory and specify video file

cd(working_dir) 

vid_file =  '20ul_test_001.avi';

disp(['Processing File: ', vid_file])
vid=VideoReader(vid_file);

%% get metaData for video processing: 

vid_dir = dir(vid_file);
vid_fps = vid.FrameRate;
num_frames = vid.NumFrames;

% gets video start and end date to construct time series:

vid_end_time = datetime([vid_dir.date]);

% create a time series of date times starting from (video end time - video duration) and ending at the video end time:

time_series = datetime(vid_end_time - (0:num_frames) * seconds(1/vid_fps), 'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
time_series = flip(time_series);  % Flips time series to be chronologically ordered
vid_start_time = datetime(time_series(1));

%% "frame by frame method" 
% *obtains SWE for each .avi file*

% preallocate variables saved in loop for speed:

h_data_cells = cell(num_frames,1);
plate_temp = nan(num_frames,1);
sum_h_area_times_dt = nan(num_frames,1);

% enter loop to process images: 

for frame_ii = 1:num_frames     
    
    frame = read(vid, frame_ii);  
    frame_gray = im2gray(frame); % convert frame of interest to gray scale
    frame_gray_cropped_wKapton = imcrop(frame_gray, colorbar_image_indexes);% crop out colorbar
    plate_temp(frame_ii) = max(max(double(frame_gray_cropped_wKapton))); % this assumes max temperature in image is the plate temperature with Kapton tape 
    frame_cropped = imcrop(frame_gray, colorbar_kapton_image_indexes); % back to orginal grayscale image... now remove colorbar and kapton tape from image
    frame_filtered = frame_cropped > min_thres; % removed below min threshold, on rbg ([0, 255]) scale 
    frame_filled = imfill(frame_filtered, 'Holes'); % clean up Hydrometeors
    frame_final = bwareaopen(frame_filled, minimum_hydro_area); % any hydrometeor whose area is less than minimum_hydro_area (set to 2 pixels) is disgarded

    % now continue on to get hydrometeor properties: 

    h_geo_prop = regionprops(frame_final, 'Area', 'Centroid'); % returns just area and centroid of each hydrometeor 
    
    % if no properties are found, go to next frame: 

    if (isempty(h_geo_prop))
        continue;
    end
    
    % build hydrometeor property matrices from regionprops values: 

    h_areaPix = cat(1, h_geo_prop.Area); % concat all hydrometoer areas in pixels
    h_centroid = round(cat(1, h_geo_prop.Centroid)); % concat all values centroid indexes
   
    % convert to length scales: 

    h_area = h_areaPix .* pix_to_m2_conversion; 
    
    
    % difference in temperature of each centroid and the plate:

    h_centroid_i = sub2ind(size(frame_cropped), h_centroid(:, 2), h_centroid(:, 1)); % find the linear index of the centriods in orginal image
    h_centroid_values = double(frame_cropped(h_centroid_i)); % intensities of centroid pixels of snow
    plate_h_dtemp = colorbar_max_temp - (h_centroid_values .* (colorbar_max_temp / plate_temp(frame_ii)));
    
    % product of hydrometeor area with the temp difference:

    h_area_times_dtemp = h_area .* plate_h_dtemp;         
    
    % sum the product of individual area and temp. diff in each frame:
    % **this is how we obtain hydrometeor mass using fbf method**

    sum_h_area_times_dt(frame_ii)=sum(h_area_times_dtemp);         
    
    % build large matrix of Hydrometeor data:

    h_data_cells{frame_ii} = cat(2, h_centroid, plate_h_dtemp, h_area); 
end

% % frame by frame mas calculation:
% 
% sum_h_area_times_dt(isnan(sum_h_area_times_dt)) = 0; % turn all NaN to 0's
% hp_area = (size(frame_cropped,1) * size(frame_cropped,2)) * pix_to_m2_conversion; % hotplate area     
% h_mass_fbf = (k_dLv*sum_h_area_times_dt) / vid_fps; % total mass evaporates in each frame 
% 
% % frame by frame SWE calculation:
% 
% SWE_fbf = h_mass_fbf / hp_area;

%% call sortPositions_v2.m 

% places snowflakes in the same row across multiple frames, making tracking possible over time

% determine the number of columns in h_data_cells to help pad matrix
% with zeros:

num_cols = [];

for ii = 1:num_frames
    if ~isempty(h_data_cells{ii})
        num_cols = size(h_data_cells{ii}, 2); % whichever cell is first detected to not be blank, grab that number of columns
        break
    end
end

% now replace empty frames with 0xnum_cols matrices:

for ii = 1:num_frames
    if isempty(h_data_cells{ii})
        h_data_cells{ii} = zeros(0, num_cols);
    end
end

% create new cell array for sorted data:

h_data_sorted = cell(size(h_data_cells));
h_data_sorted{1} = h_data_cells{1};

% loop over every data cell corresponding to each frame:
 
for frame_jj = 2:num_frames
    h_data_sorted{frame_jj} = h_data_cells{frame_jj};
    h_data_sorted{frame_jj} = sortPositions_v2(h_data_sorted{frame_jj-1}, h_data_sorted{frame_jj}, sort_threshold);
end


% sum all properties for frames where more than one hydrometeor was found:
% (this is when a droplet splits into multiple particles)

h_data_sorted = cellfun(@(x) sum(x,1), h_data_sorted, 'UniformOutput', false);

% return frame with max number of hydrometeors:

max_h_obs = max(cellfun(@(x) size(x, 1), h_data_sorted, 'UniformOutput', 1));

% now pad the data with zeros so the frames all have the same number:

h_data_sorted = cellfun(@(x) cat(1, x, zeros(max_h_obs - size(x, 1), width(h_data_sorted{1}))),...
    h_data_sorted, 'UniformOutput', 0);

%% isolating the variables and put them into a matrix to work with

% for reference: [h_centroid(1), h_centroid(2), plate_h_dtemp, h_area]

dT_fbf = cellfun(@(x) x(:, 3), h_data_sorted, 'UniformOutput', 0);
area_fbf = cellfun(@(x) x(:, 4), h_data_sorted, 'UniformOutput', 0);

% convert to matrix: 

dT_fbf = cat(2,dT_fbf{:});
area_fbf = cat(2,area_fbf{:}); % snowflake area 

%% "particle by particle method" 

% this method uses the data obtained from fbf, but now isolates each
% snowflake per frame 

% initialize parameters:

all_h_appears_ind = []; % initial time index of each snowflake in time array
h_delta_temp = {}; % hydrometeor delta temp over time
h_area = {}; % hydrometeor area over time

h_max_area = {}; % hydrometeor maximum area
h_delta_time = {}; % hydrometeor time to melt and evaporate
h_delta_temp_max = {}; % hydrometeor centroid's max intensity
h_delta_temp_mean = {}; % hydrometeor centroid's mean intensity
snowflakeAreaTimesDeltaTemp = {}; % summing area*deltaTemp for each snowflake

% loop over all Hydrometeors for each frame: 

for h_ii = 1:max_h_obs

    h_appear_evap_bool = diff(area_fbf(h_ii, :) > 0); % creates an array 
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

    h_dT_tmp = dT_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
    h_area_tmp = area_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
    
    % isolate for positive values:

    h_area_tmp_bool = h_area_tmp > 0; 
    
    % finds how many 'continuous' snowflakes there are: indices of contiguous '1s':

    propstemp = regionprops(h_area_tmp_bool, 'PixelIdxList'); 

    % grab properties of each snowflake over their 'life': 

        for jj = 1:numel(propstemp)

            h_delta_temp{end+1} = h_dT_tmp(propstemp(jj).PixelIdxList); % hydrometeor delta temp over time
            h_area{end+1} = h_area_tmp(propstemp(jj).PixelIdxList); % hydrometeor area over time
            h_max_area{end+1} = max(h_area_tmp(propstemp(jj).PixelIdxList)); % max snowflake area over time
            h_delta_time{end+1} = numel(propstemp(jj).PixelIdxList); % time of snowflake's life 
            h_delta_temp_max{end+1} = max(h_dT_tmp(propstemp(jj).PixelIdxList)); % max temperature difference between snowflake and plate 
            h_delta_temp_mean{end+1} = mean(h_dT_tmp(propstemp(jj).PixelIdxList)); % mean temperature difference between snowflake and plate
                            
            % multiply the area by the temperature difference for that hydrometeor:

            h_area_times_dT_pbp = h_area_tmp(propstemp(jj).PixelIdxList) .* h_dT_tmp(propstemp(jj).PixelIdxList);
            snowflakeAreaTimesDeltaTemp{end+1} = sum(h_area_times_dT_pbp); % sum over the snowflakes life and store 

        end
    end
end

%% convert to matrixes

h_delta_temp = cell2mat(h_delta_temp);
h_area = cell2mat(h_area);
snowflakeAreaTimesDeltaTemp = cell2mat(snowflakeAreaTimesDeltaTemp); % used to calculate mass (calibrate k/d/L_v)
h_max_area=cell2mat(h_max_area); % actual area
h_delta_temp_max=cell2mat(h_delta_temp_max); % temperature diffrence between plate and water droplet using max intensity
h_delta_temp_mean=cell2mat(h_delta_temp_mean); % temperature diffrence between plate and water droplet using mean intensity 

%% conversions and calculations using PBP data

% solve for k_dLv coefficient assuming that each droplet was 20uL:

k_dLv_calib = (mass_calib * vid_fps) / median(snowflakeAreaTimesDeltaTemp); 

fprintf('Calibrated k_dLv = %.6e\n', k_dLv_calib);

% now use that coefficent value to calculate mass and volume of droplet

m_pbp = (k_dLv_calib .* snowflakeAreaTimesDeltaTemp) / vid_fps; % kg
v_pbp = 1e6 * m_pbp; % uL

% deviation from nominal 20 uL

deltaV = v_pbp - 20;
pct_off  = 100 * deltaV / 20;  

T = table( ...
    m_pbp(:), ...
    v_pbp(:), ...
    deltaV(:), ...
    pct_off(:), ...
    'VariableNames', { ...
        'Mass_kg', ...
        'Volume_uL', ...
        'DeltaVolume_uL', ...
        'PercentOff_20uL' ...
    } ...
);
