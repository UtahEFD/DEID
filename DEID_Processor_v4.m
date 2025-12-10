%% DEID AVI File Processing Code

% Outputs filtered, unfiltered, and time-averaged particle-by-particle .csv
% files for many .avi files using a parallelized for loop. For each .avi file,
% a summary file is also generated. 

% AUTHOR : Dhiraj Singh, Benjamin Silberman, Travis Morrison, Alex Blackmer
                                       
clear, clc, close all
delete(gcp('nocreate')); 
%% set filepath, output directory, and file name for saving  

working_dir = '/uufs/chpc.utah.edu/common/home/snowflake4/DEID_files/2024_2025/apr01_storm';
output_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/testData/new';
storm_output = '_apr0125';

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

min_thres = 70; % minimum threshold number in image accepted rbg ([0 255]) scale
minimum_hydro_area = 2; % the minimum number of pixels a hydrometeor must contain to be analyzed 
colorbar_max_temp = 145; % max temperature set in colorbar on the physical screen of the tir software
sort_threshold = 20; % this is the RMS threshold between succesive images of snowflakes used in the sortPostitions_v2. dhiraj calibrated this in the lab. 
minimum_drop_life = 0; % minimum number of frames a drop has to be visable to be processed
areaTol = 0; 
SWEfactor_threshold = 1.85; % maximum value of tolerable SWE factor
evapTime_min = 1/15; % minimum time a snowflake has to appear on hotplate to be processed
evapTime_max = 30; % maximum time a snowflake can appear on hotplate to be processed 
noiseThresh = 1000; % # of times a centroid can appear in an .avi file to be considered real 

% DEID specific parameters:

colorbar_image_indexes = [1 1 384 288]; % location of colorbar in pixel locations
crop_index = 43; % use this to specify indices to crop out kapton tape
colorbar_kapton_image_indexes = [1 (colorbar_image_indexes(2)+crop_index) 383 (colorbar_image_indexes(4)-crop_index)]; % location of Kapton tape in pixel locations
k_dLv = 0.0035; % calibration constant; in paper, thermal conductivity (k) of water. See sect 4.1 in Dihiraj's paper -> (k/d(_eff))/Latent heat of vaporazation [units?]
l_constant = 2.594e06; % latent heat of vaporazation of water, should be a function of tempertaure (Look at Stull textbook) [J/kg]
% Eqn. (13) in Dhiraj's density paper: c = (L_vv) / (L_ff*C_melt)
% c = hf_rho_coeff
hf_rho_coeff = 1.01e05; % [K*s*m^-1]

%% move to working directory and identify video files 

% get all directory contents:

cd(working_dir);
directory = dir("."); 

% loop through directory items:

file_names = cell(1, length(directory));
vid_date   = cell(1, length(directory));
count = 0;
for file_i = 1:length(directory)
    if ~directory(file_i).isdir && endsWith(directory(file_i).name, '.avi', 'IgnoreCase', true) % check if it's a file and ends with .avi
        count = count + 1;
        file_names{count} = directory(file_i).name; % store filename
        vid_dir = dir(directory(file_i).name);
        vid_date{count} = datetime(vid_dir.date);  % get timestamp from the file metadata
    end
end

% trim unused cells: 

file_names = file_names(1:count);
vid_date   = vid_date(1:count);

% sort filenames by date:

[~, sort_idx] = sort([vid_date{:}]);
file_names = file_names(sort_idx);

% **when testing**
% file_names = {'Mar06_2023_29.avi'};

%% begin DEID video processing
% Parfor loop parallelizes processing by distributing each video file to a
% Matlab worker on each CPU core. 

% preallocate variables:

pbp_table_cell = cell(length(file_names),1);
pbp_table_filtered_cell = cell(length(file_names),1); 
avi_summary_table_cell = cell(length(file_names),1);
pbp_table_retimed_cell = cell(length(file_names),1); 

parpool(10); 
parfor file_i = 1:length(file_names)

    filename = file_names{file_i};
    disp(['Processing File: ', filename])
    vid=VideoReader(filename);
    
    % get metaData for video processing: 

    vid_dir = dir(filename);
    vid_fps = vid.FrameRate;
    num_frames = vid.NumFrames;

    % gets video start and end date to construct time series:

    vid_end_time = datetime([vid_dir.date]);
    
    % create a time series of date times starting from (video end time - video duration) and ending at the video end time:

    time_series = datetime(vid_end_time - (0:num_frames) * seconds(1/vid_fps), 'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
    time_series = flip(time_series);  % flips time series to be chronologically ordered
    vid_start_time = datetime(time_series(1));

    %% identify noisy hydrometeors using centroids:
    
    % store full video in CHPC RAM:

    frames = cell(num_frames,1);
    
    for ii = 1:num_frames
        frames{ii} = read(vid,ii);
    end
    
    % collect all centroids:

    allCentroids = [];
    
    for ii = 1:num_frames
        frame_gray = im2gray(frames{ii});
        frame_cropped = imcrop(frame_gray, colorbar_kapton_image_indexes);
        frame_filtered = frame_cropped > min_thres;
        frame_filled   = imfill(frame_filtered,'holes');
        frame_cleaned  = bwareaopen(frame_filled, minimum_hydro_area);
    
        allProps = regionprops(frame_cleaned,'Centroid');
        if ~isempty(allProps)
            allCentroids = [allCentroids; cat(1,allProps.Centroid)];
        end

    end

    % if the video file is blank, skip entirely:
    % 
    % if isempty(allCentroids)
    % 
    %     warning('Video %s contains no detectable hydrometeors. Skipping.', filename);
    % 
    %     pbp_table = timetable();
    %     pbp_table_filtered = timetable();
    %     avi_summary_table = timetable();
    %     pbp_table_retimed = timetable();
    % 
    %     pbp_table_cell{file_i} = pbp_table;
    %     pbp_table_filtered_cell{file_i} = pbp_table_filtered;
    %     avi_summary_table_cell{file_i} = avi_summary_table;
    %     pbp_table_retimed_cell{file_i} = pbp_table_retimed;
    % 
    %     continue
    % 
    % end
    
    % identify noisy centroids:

    [uniqueC,~,idxC] = unique(allCentroids,'rows');
    counts = accumarray(idxC,1);
    noiseMask   = counts > noiseThresh;
    noiseCentroids = uniqueC(noiseMask,:);
    
    %% optional: view ten most repeated centroids:
    % [counts_sorted, sortIdx] = sort(counts, 'descend');
    % top_centroids = uniqueC(sortIdx(1:15), :);
    % top_counts = counts_sorted(1:15);
    % table([top_centroids, top_counts])

    %% "frame by frame method"; this is how we obtain SWE for each .avi file
    
    % preallocate variables saved in loop for speed:

    h_data_cells = cell(num_frames,1);
    plate_temp = nan(num_frames,1);
    noisyA = cell(num_frames,1); 
    sum_h_area_times_dt = nan(num_frames,1);
    
    % enter loop to process images: 
    
    for frame_ii = 1:num_frames     
        frame = frames{frame_ii};  
        frame_gray = im2gray(frame); % convert frame of interest to gray scale
        frame_gray_cropped_wKapton = imcrop(frame_gray, colorbar_image_indexes);% crop out colorbar
        plate_temp(frame_ii) = max(max(double(frame_gray_cropped_wKapton))); % this assumes max temperature in image is the plate temperature with Kapton tape 
        frame_cropped = imcrop(frame_gray, colorbar_kapton_image_indexes); % back to orginal grayscale image... now remove colorbar and kapton tape from image
        frame_filtered = frame_cropped > min_thres; % removed below min threshold, on rbg ([0, 255]) scale 
        frame_filled = imfill(frame_filtered, 'Holes'); % clean up Hydrometeors
        frame_final = bwareaopen(frame_filled, minimum_hydro_area); % any hydrometeor whose area is less than minimum_hydro_area (set to 2 pixels) is disgarded
        
        % remove centroids that appear more than num_frames/4 times:

        props = regionprops(frame_final, 'Area', 'Centroid','PixelIdxList');

        if ~isempty(props)

            frame_centroids = cat(1, props.Centroid); % collect centroids of particles on frame 
            frame_area = cat(1, props.Area); % collect areas of particles on frame 

        if isempty(noiseCentroids)

            isNoise = false(size(frame_centroids,1),1); % basically do nothing

        else

            isNoise = ismember(frame_centroids, noiseCentroids, 'rows'); % logical array identifying which centroids on the frame are noisy ones 

        end

            noisyA{frame_ii} = sum(frame_area(isNoise, :)); % store area to subtract from hpArea later 

            % black out noisy hydrometeors:

            for k = find(isNoise)'
                frame_final(props(k).PixelIdxList) = 0;
            end
    
        end

        % now continue on to get hydrometeor properties: 

        h_geo_prop = regionprops(frame_final, 'PixelIdxList', 'MajorAxisLength', 'MinorAxisLength', 'Centroid', 'Area', 'BoundingBox', 'Perimeter'); % returns the centroid, the area , and the bounding box (left, top, width, height) of each blob
        
        % if no properties are found, go to next frame: 

        if (isempty(h_geo_prop))
            continue;
        end

        % PCA-BASED CIRCUMSCRIBED ELLIPSE AREA PER HYDROMETEOR
        
        h_PCAellipseAreaM = zeros(length(h_geo_prop),1);
        
        for ii = 1:length(h_geo_prop)
        
            % extract pixel coordinates of hydrometeor:
            
            pixList = h_geo_prop(ii).PixelIdxList;
            [r, c] = ind2sub(size(frame_final), pixList);
            pts = [c, r];  % nx2 array of pixel coordinates for each hydrometeor 
        
            % center the points:
            C    = mean(pts,1);   % centroid of hydrometeor in pixel coordinates 
            pts0 = pts - C;       % centered coordinates of pixels 
        
            % covariance + eigen decomposition:
            Sigma    = cov(double(pts0));
            [V, D]   = eig(Sigma);

            % sort eigenvectors to ensure major & minor axes:
            [~, idx] = sort(diag(D), 'descend');
            V = V(:, idx);
            D = diag(sort(diag(D), 'descend'));
        
            % rotate points into PCA basis:
            ptsRot = pts0 * V;
        
            % first: maximum absolute extent in PCA axes (inscribed semi-axes):
            a0_pix = max(abs(ptsRot(:,1)));   % initial semi-major axis
            b0_pix = max(abs(ptsRot(:,2)));   % initial semi-minor axis
        
            % scale ellipse so it CIRCUMSCRIBES all points:
            normVals = (ptsRot(:,1)/a0_pix).^2 + (ptsRot(:,2)/b0_pix).^2;
            s        = sqrt(max(normVals));   % >= 1
        
            % final circumscribing semi-axes (keep same variable names):
            a_pix = s * a0_pix + 0.5;
            b_pix = s * b0_pix + 0.5;
        
            % area of ellipse (pixel units) using circumscribing semi-axes:
            ellipse_area_pix = pi * a_pix * b_pix;
        
            % convert to m^2:
            h_PCAellipseAreaM(ii) = ellipse_area_pix * pix_to_m2_conversion;

            % checking complexity:
            areaPix = numel(pixList); 
            Cx_pix = ellipse_area_pix / areaPix;
            if Cx_pix < 1
                fprintf('Frame %d, hydro %d: Cx_pix = %.3f\n', frame_ii, ii, Cx_pix);
            end
        
        end

        
        % build hydrometeor property matrices from regionprops values: 

        h_bounding_box = cat(1,h_geo_prop.BoundingBox); % concat all values to bounding box indices in pixels
        rect_widthPix = h_bounding_box(:,3); 
        rect_heightPix = h_bounding_box(:,4);
        h_centroid = round(cat(1, h_geo_prop.Centroid)); % concat all values centroid indexes
        h_perimeterPix = cat(1,h_geo_prop.Perimeter); % concat all values hydrometeor perimeters in pixels
        h_areaPix = cat(1, h_geo_prop.Area); % concat all hydrometoer areas in pixels
        h_majorPix = cat(1,h_geo_prop.MajorAxisLength); % concat all major axis in pixels
        h_minorPix = cat(1,h_geo_prop.MinorAxisLength); % concat all minor axis in pixels
        
        % convert to length scales: 

        rect_widthM = rect_widthPix * pix_to_m_conversion;
        rect_heightM = rect_heightPix * pix_to_m_conversion;
        h_perimeterM = h_perimeterPix * pix_to_m_conversion; 
        h_area = h_areaPix .* pix_to_m2_conversion; 
        h_majorM = h_majorPix * pix_to_m_conversion;
        h_minorM = h_minorPix * pix_to_m_conversion;

        % calculate circumscribed areas: 

        h_rectAreaM = rect_widthM .* rect_heightM; % rectangle
        h_circleAreaM = (pi * h_majorM.^2)/4; % circumscribed circle using major axis 
        
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

        h_data_cells{frame_ii} = cat(2, h_centroid, plate_h_dtemp, h_perimeterM, h_area, h_rectAreaM, h_circleAreaM, h_PCAellipseAreaM, rect_widthM, rect_heightM, h_majorM); 
    end

    % frame by frame SWE calculation:

    sum_h_area_times_dt(isnan(sum_h_area_times_dt)) =0; % turn all NaN to 0's
    hp_area = ((size(frame_cropped,1) * size(frame_cropped,2)) - mean([noisyA{:}])) * pix_to_m2_conversion; % hotplate area     
    h_mass_fbf = (k_dLv*sum_h_area_times_dt) / vid_fps; % total mass evaporates in each frame
    h_mass_fbf_min = min(h_mass_fbf); % we know the plate should be empty when it is not snowing..
    h_mass_fbf = h_mass_fbf - h_mass_fbf_min; % subtract off min mass on a frame to account for any resiude
    SWE_fbf = h_mass_fbf / hp_area;

    % find the minimum SWE in all frames within a video, and subtract from
    % SWE (way of handling residue):

    time_series_fbf = time_series(1:length(SWE_fbf));
    
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

    % but if all frames are empty, skip: 
    if isempty(num_cols)
        warning('Video %s contains no hydrometeors. Skipping.', filename);
    
        pbp_table = timetable();
        pbp_table_filtered = timetable();
        avi_summary_table = timetable();
        pbp_table_retimed = timetable();
    
        pbp_table_cell{file_i} = pbp_table;
        pbp_table_filtered_cell{file_i} = pbp_table_filtered;
        avi_summary_table_cell{file_i} = avi_summary_table;
        pbp_table_retimed_cell{file_i} = pbp_table_retimed;
    
        continue
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

    % return frame with max number of hydrometeors:

    max_h_obs = max(cellfun(@(x) size(x, 1), h_data_sorted, 'UniformOutput', 1));

    % now pad the data with zeros so the frames all have the same number:

    h_data_sorted = cellfun(@(x) cat(1, x, zeros(max_h_obs - size(x, 1), width(h_data_sorted{1}))),...
        h_data_sorted, 'UniformOutput', 0);
    
    %% isolating the variables and put them into a matrix to work with
    
    % for reference: [h_centroid(1), h_centroid(2), plate_h_dtemp,... 
    % h_perimeterM, h_area, h_rectAreaM, h_circleAreaM,... 
    % h_PCAellipseAreaM, rect_widthM, rect_heightM, h_majorM]
    
    dT_fbf = cellfun(@(x) x(:, 3), h_data_sorted, 'UniformOutput', 0);
    perimeter_fbf = cellfun(@(x) x(:, 4), h_data_sorted, 'UniformOutput', 0);
    area_fbf = cellfun(@(x) x(:, 5), h_data_sorted, 'UniformOutput', 0);
    rectArea_fbf = cellfun(@(x) x(:, 6), h_data_sorted, 'UniformOutput', 0);
    circleArea_fbf = cellfun(@(x) x(:, 7), h_data_sorted, 'UniformOutput', 0);
    ellipseArea_fbf = cellfun(@(x) x(:, 8), h_data_sorted, 'UniformOutput', 0); 
    rectWidth_fbf = cellfun(@(x) x(:, 9), h_data_sorted, 'UniformOutput', 0);
    rectHeight_fbf = cellfun(@(x) x(:, 10), h_data_sorted, 'UniformOutput', 0);
    h_majorAxis_fbf = cellfun(@(x) x(:, 11), h_data_sorted, 'UniformOutput', 0);
   
    % convert to matrix: 
    
    dT_fbf = cat(2,dT_fbf{:});
    perimeter_fbf = cat(2, perimeter_fbf{:}); % snowflake perimeter 
    area_fbf = cat(2,area_fbf{:}); % snowflake area 
    rectArea_fbf = cat(2,rectArea_fbf{:}); % circumscribed rectangle area 
    circleArea_fbf = cat(2,circleArea_fbf{:}); % circumscribed circle area using majorAxis as D 
    ellipseArea_fbf = cat(2, ellipseArea_fbf{:}); % circumscribed ellipse area using PCA 
    ellipseArea_fbf(isnan(ellipseArea_fbf)) = 0; 
    rectWidth_fbf = cat(2,rectWidth_fbf{:}); % circumscribed rectangle width
    rectHeight_fbf = cat(2,rectHeight_fbf{:}); % circumscribed rectangle height 
    h_majorAxis_fbf = cat(2, h_majorAxis_fbf{:}); % hydrometeor major axis 

    %% "particle by particle method" 

    % this method uses the data obtained from fbf, but now isolates each
    % snowflake per frame 
    
    % initialize parameters:

    all_h_appears_ind = []; % initial time index of each snowflake in time array
    h_delta_temp = {}; % hydrometeor delta temp over time
    h_perimeter = {}; % hydrometeor perimeter over time
    h_area = {}; % hydrometeor area over time
    h_rectArea = {}; % circumscribed rectangle area
    h_circleArea = {}; % circumscribed circle area
    h_ellipseArea = {}; % circumscribed ellipse area 
    h_rectWidth = {}; % circumscribed rectangle width
    h_rectHeight = {}; % circumscribed rectangle height
    h_majorAxis = {}; % hydrometeor max axis length  
    
    deltaTemp_range = []; % difference between max and min values of hydrometeor deltaTemp
    deltaTemp_residue_flags = []; % flags for when deltaTemp does not change
    hArea_range = []; % difference between max and min values of hydrometeor area
    hArea_residue_flags = []; % flags for when area does not change 
    
    h_max_area = {}; % hydrometeor maximum area
    h_delta_time = {}; % hydrometeor time to melt and evaporate
    h_delta_temp_max = {}; % hydrometeor centroid's max intensity
    h_delta_temp_mean = {}; % hydrometeor centroid's mean intensity
    h_mass_pbp = {}; % hydrometeor calculated mass

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
        h_perimeter_tmp = perimeter_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
        h_area_tmp = area_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
        h_rectArea_tmp = rectArea_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
        h_circleArea_tmp = circleArea_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1); 
        h_ellipseArea_tmp = ellipseArea_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1); 
        h_rectWidth_tmp = rectWidth_fbf(h_ii,h_appears_ind(1):h_evaps_ind(end)+1);
        h_rectHeight_tmp = rectHeight_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
        h_majorAxis_tmp = h_majorAxis_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
        
        % isolate for positive values:

        h_area_tmp_bool = h_area_tmp > 0; 
        
        % finds how many 'continuous' snowflakes there are: indices of contiguous '1s':

        propstemp = regionprops(h_area_tmp_bool, 'PixelIdxList'); 

        % grab properties of each snowflake over their 'life': 

            for jj = 1:numel(propstemp)
                h_delta_temp{end+1} = h_dT_tmp(propstemp(jj).PixelIdxList); % hydrometeor delta temp over time
                h_area{end+1} = h_area_tmp(propstemp(jj).PixelIdxList); % hydrometeor area over time
                
                h_perimeter{end+1} = max(h_perimeter_tmp(propstemp(jj).PixelIdxList)); % max hydrometeor perimeter over time
                h_rectArea{end+1} = max(h_rectArea_tmp(propstemp(jj).PixelIdxList)); % max rectangle area over time
                h_circleArea{end+1} = max(h_circleArea_tmp(propstemp(jj).PixelIdxList)); % max circle area over time
                h_ellipseArea{end+1} = max(h_ellipseArea_tmp(propstemp(jj).PixelIdxList)); % max ellipse area over time 
                h_rectWidth{end+1} = max(h_rectWidth_tmp(propstemp(jj).PixelIdxList)); % max rectangle width over time
                h_rectHeight{end+1} = max(h_rectHeight_tmp(propstemp(jj).PixelIdxList)); % max rectangle height over time
                h_majorAxis{end+1} = max(h_majorAxis_tmp(propstemp(jj).PixelIdxList)); % max major axis over time
                h_max_area{end+1} = max(h_area_tmp(propstemp(jj).PixelIdxList)); % max snowflake area over time
                h_delta_time{end+1} = numel(propstemp(jj).PixelIdxList); % time of snowflake's life 
                h_delta_temp_max{end+1} = max(h_dT_tmp(propstemp(jj).PixelIdxList)); % max temperature difference between snowflake and plate 
                h_delta_temp_mean{end+1} = mean(h_dT_tmp(propstemp(jj).PixelIdxList)); % mean temperature difference between snowflake and plate
                
                % filtering:

                deltaTemp_range(end+1) = range(h_delta_temp{end});
                deltaTemp_residue_flags(end+1) = (range(h_delta_temp{end}) == 0); % flag snowflakes whose delta temp does not change over time
                hArea_range(end+1) = range(h_area{end});
                hArea_residue_flags(end+1) = (range(h_area{end}) <= areaTol); % flag snowflakes whose area does not change over time
                                
                % multiply the area by the temperature difference for that hydrometeor:

                h_area_times_dT_pbp = h_area_tmp(propstemp(jj).PixelIdxList) .* h_dT_tmp(propstemp(jj).PixelIdxList);

                % now integrate (sum over snowflake's life):

                h_mass_pbp{end+1} = sum(k_dLv * h_area_times_dT_pbp); % total mass that each snoflake contains across its entire life cylce

            end
        end
    end

    %% convert to matrixes
    
    h_delta_temp = cell2mat(h_delta_temp);
    h_perimeter = cell2mat(h_perimeter); 
    h_area = cell2mat(h_area);
    h_rectArea=cell2mat(h_rectArea);
    h_circleArea = cell2mat(h_circleArea);
    h_ellipseArea = cell2mat(h_ellipseArea); 
    h_rectWidth=cell2mat(h_rectWidth); 
    h_rectHeight=cell2mat(h_rectHeight); 
    h_majorAxis = cell2mat(h_majorAxis);

    h_mass_pbp=cell2mat(h_mass_pbp); % hydrometeor mass
    h_max_area=cell2mat(h_max_area); % actual area
    h_delta_temp_max=cell2mat(h_delta_temp_max); % temperature diffrence between plate and water droplet using max intensity
    h_delta_temp_mean=cell2mat(h_delta_temp_mean); % temperature diffrence between plate and water droplet using mean intensity 
    
    %% conversions and calculations using PBP data

    h_mass_pbp = h_mass_pbp / vid_fps; 
    h_dEff = ((4/pi) * h_max_area).^(1/2); % convert area to diameter of hydrometeor
    h_circlePerimeter = pi*h_dEff; % circumscribed circle perimeter around a hydrometeor 
    h_evap_time = cell2mat(h_delta_time) * (1 / vid_fps); % evaporation time
    h_vol_sph = (3/4) * h_max_area.^(3/2); % spherical volume
    h_rho_sph = h_mass_pbp ./ h_vol_sph; % density calculation: spherical assumption
    h_energy_per_time = h_mass_pbp * l_constant ./ (h_max_area .* h_evap_time); % heat flux method: energy per unit area per time
    h_height = h_evap_time .* h_delta_temp_mean; 
    h_rho_hfd = (hf_rho_coeff * h_mass_pbp) ./ (h_max_area .* h_evap_time .* h_delta_temp_mean); 
    h_vol_hfd = h_mass_pbp ./ h_rho_hfd; % volume of each snowflakes using mean heat flux method density
    h_initial_time_indices = all_h_appears_ind(1:length(h_mass_pbp));
    h_initial_time = time_series(h_initial_time_indices);
    eqWaterDrop_area = ((9*pi)/16)^(1/3) * (h_mass_pbp ./ rho_water).^(2/3); % surface area of equivalent water droplet (See POF Singh et al. 2023)
    
    % Cx and SDI calculations:

    complexity1 = h_rectArea ./ h_max_area; % complexity using boundingBox area (See CRST Morrison et al. 2023)
    complexity2 = h_circleArea ./ h_max_area; % complexity using circle area ((pi*majorAxis^2)/4)
    complexity3 = h_circlePerimeter ./ h_perimeter; % complexity using circle perimeter and perimeter of snowflake
    complexity4 = h_ellipseArea ./ h_max_area; % complexity using ellipse area (PCA) 
    
    sdi = h_max_area ./ eqWaterDrop_area ; % SDI (See CRST Morrison et al. 2023)
    
    % SWE and snow calculations:

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

    % store PBP data as table for post processing:

    pbp_table = table(h_initial_time', h_evap_time', ... 
        h_mass_pbp', h_dEff',  h_perimeter', h_max_area', h_rectArea', ... 
        h_circleArea', h_ellipseArea', eqWaterDrop_area', h_rho_sph', h_rho_hfd', h_vol_hfd', h_vol_sph', ... 
        h_delta_temp_max', h_delta_temp_mean', complexity1', complexity2', complexity3', complexity4', ... 
        sdi', SWE_pbp', SWE_fbf_particles', SWE_pbp_accumulated', ... 
        SWE_fbf_accumulated', snow_pbp', snow_fbf', snow_pbp_accumulated', ... 
        snow_fbf_accumulated', SWE_factor_particles', deltaTemp_range', hArea_range', deltaTemp_residue_flags', hArea_residue_flags');
    pbp_table.Properties.VariableNames = {'Time', 'Evap Time (s)', 'Mass (kg)', 'Eff Diameter (m)',  ...
        'Perimeter (m)', 'Snowflake Area (m^2)', 'Rectangle Area (m^2)', 'Circle Area (m^2)',... 
        'Ellipse Area (m^2)', 'Water Droplet Area (m^2)', 'Spherical Density (kg/m^3)', 'Heat Flux Density (kg/m^3)', ... 
        'Heat Flux Volume (m^3)', 'Spherical Volume (m^3)', 'Delta Temp Max', 'Delta Temp Mean', ...
        'Rect Complexity', 'Circle Complexity', 'Perimeter Complexity', 'Ellipse Complexity', 'SDI','PBP SWE (mm)','FBF SWE (mm)', 'PBP SWE Accumulation (mm)', ...
        'FBF SWE Accumulation (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'PBP Snow Accumulation (mm)', ...
        'FBF Snow Accumulation (mm)', 'SWE factor', 'Delta Temp Range', 'Area Range', 'Delta Temp Flag', 'Area Flag'};
    pbp_table = table2timetable(pbp_table); % convert to timetable 
    pbp_table = sortrows(pbp_table, 'Time'); % sort by time 
    pbp_table.("PBP SWE Accumulation (mm)") = cumsum(pbp_table.("PBP SWE (mm)"));
    pbp_table.("FBF SWE Accumulation (mm)") = cumsum(pbp_table.("FBF SWE (mm)"));
    pbp_table.("PBP Snow Accumulation (mm)") = cumsum(pbp_table.("PBP Snow (mm)"));
    pbp_table.("FBF Snow Accumulation (mm)") = cumsum(pbp_table.("FBF Snow (mm)")); 
    pbp_table.("Missing Data") = false(height(pbp_table),1); % add a flag column to distinguish missing .avi data

    %% filtering starts here

    % filters data for any hydrometeor whose area doesn't change over it's
    % life cycle, and whose deltaTemp doesn't change over it's life cycle:

    pbp_table_filtered = pbp_table(pbp_table.('Area Flag') ~= 1 &... 
        pbp_table.('Delta Temp Flag') ~= 1 &...
        pbp_table.("Evap Time (s)") > evapTime_min &...
        pbp_table.("Evap Time (s)") < evapTime_max, :);

    % re-accumulate totals: 

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
        perRow = mean(prev_data.('Perimeter (m)'));
        areaRow =  mean(prev_data.('Snowflake Area (m^2)'));
        rectAreaRow = mean(prev_data.("Rectangle Area (m^2)")); 
        circAreaRow = mean(prev_data.("Circle Area (m^2)"));
        ellipseAreaRow = mean(prev_data.("Ellipse Area (m^2)"));
        waterDropletAreaRow = mean(prev_data.("Water Droplet Area (m^2)")); 
        volumeHFDrow = sum(prev_data.("Heat Flux Volume (m^3)"));
        volumeSPHrow = sum(prev_data.("Spherical Volume (m^3)"));
        deltaTempMaxRow = mean(prev_data.("Delta Temp Max"));
        deltaTempMeanRow = mean(prev_data.("Delta Temp Mean"));
        cx1Row = mean(prev_data.('Rect Complexity'));
        cx2Row = mean(prev_data.('Circle Complexity'));
        cx3Row = mean(prev_data.('Perimeter Complexity'));
        cx4Row = mean(prev_data.('Ellipse Complexity'));
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
        
        % compile into a timetable:

        new_row = table(timeRow,... 
            evapTimeRow, massRow, effDiaRow, perRow, areaRow, rectAreaRow, circAreaRow, ellipseAreaRow, waterDropletAreaRow, ...
            densitySPHrow, densityHFDrow, volumeHFDrow, volumeSPHrow, deltaTempMaxRow, ...
            deltaTempMeanRow, cx1Row, cx2Row, cx3Row, cx4Row, sdiRow, PBPsweRow, FBFsweRow, ...
            sum(pbp_table.("PBP SWE (mm)"))+PBPsweRow,...
            sum(pbp_table.("FBF SWE (mm)"))+FBFsweRow, ...
            PBPsnowRow, FBFsnowRow, ...
            sum(pbp_table.("PBP Snow (mm)"))+PBPsnowRow, ...
            sum(pbp_table.("FBF Snow (mm)"))+FBFsnowRow, SWEfactorRow, ...
            tempRangeRow, aRangeRow, tempFlagRow, aFlagRow, ...
            'VariableNames', {'Time', 'Evap Time (s)', 'Mass (kg)', 'Eff Diameter (m)',  ...
        'Perimeter (m)', 'Snowflake Area (m^2)', 'Rectangle Area (m^2)', 'Circle Area (m^2)',... 
        'Ellipse Area (m^2)', 'Water Droplet Area (m^2)', 'Spherical Density (kg/m^3)', 'Heat Flux Density (kg/m^3)', ... 
        'Heat Flux Volume (m^3)', 'Spherical Volume (m^3)', 'Delta Temp Max', 'Delta Temp Mean', ...
        'Rect Complexity', 'Circle Complexity', 'Perimeter Complexity', 'Ellipse Complexity', 'SDI','PBP SWE (mm)','FBF SWE (mm)', 'PBP SWE Accumulation (mm)', ...
        'FBF SWE Accumulation (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'PBP Snow Accumulation (mm)', ...
        'FBF Snow Accumulation (mm)', 'SWE factor', 'Delta Temp Range', 'Area Range', 'Delta Temp Flag', 'Area Flag'});
        
        % assign a time and logical value for missing data to the new row:
        
        new_row.('Missing Data') = true;
        
        % convert to timetable:
        
        new_row = table2timetable(new_row); 
        
        % now append new row to particle_output_table:
        
        pbp_table = [pbp_table; new_row];

        %%  create a summary table with actual data 
        
        avi_summary_table = table(pbp_table_filtered.Time(1));
        avi_summary_table.duration = (pbp_table_filtered.Time(end)-pbp_table_filtered.Time(1)); 
        avi_summary_table.rectCx = mean(pbp_table_filtered.('Rect Complexity'));
        avi_summary_table.circleCx = mean(pbp_table_filtered.('Circle Complexity'));
        avi_summary_table.perCx = mean(pbp_table_filtered.('Perimeter Complexity'));
        avi_summary_table.ellipCx = mean(pbp_table_filtered.('Ellipse Complexity')); 
        avi_summary_table.sdi = mean(pbp_table_filtered.SDI);
        avi_summary_table.rho = sum(pbp_table_filtered.("Mass (kg)")) / sum(pbp_table_filtered.("Heat Flux Volume (m^3)"));
        avi_summary_table.pbpSWE = pbp_table_filtered.("PBP SWE Accumulation (mm)")(end);
        avi_summary_table.fbfSWE = pbp_table_filtered.("FBF SWE Accumulation (mm)")(end);
        avi_summary_table.pbpSnow = pbp_table_filtered.("PBP Snow Accumulation (mm)")(end);
        avi_summary_table.fbfSnow = pbp_table_filtered.("FBF Snow Accumulation (mm)")(end);
        avi_summary_table.hotPlateArea = hp_area; 
        avi_summary_table.SWEfactor = SWE_factor; 
        avi_summary_table.minMass = h_mass_fbf_min;

        %% create a table with 10 minute averaged data

        % data to average:

        avg_cols = {'Rect Complexity', 'Circle Complexity', 'Perimeter Complexity', 'Ellipse Complexity', 'SDI', 'Eff Diameter (m)', 'Snowflake Area (m^2)', 'Evap Time (s)', 'SWE factor'};
        avg_table = retime(pbp_table_filtered(:, avg_cols), 'regular', 'mean', 'TimeStep', time_step);

        % data to sum:

        sum_cols = {'Mass (kg)', 'Heat Flux Volume (m^3)'};
        sum_table = retime(pbp_table_filtered(:, sum_cols), 'regular', 'sum', 'TimeStep', time_step);
        
        % join tables:

        pbp_table_retimed = horzcat(avg_table, sum_table);
        
        % calculate density:
        
        pbp_table_retimed.('Heat Flux Density (kg/m^3)') = pbp_table_retimed.('Mass (kg)')./ pbp_table_retimed.('Heat Flux Volume (m^3)');
        
        % calculate SWE:
        
        pbp_table_retimed.('FBF SWE (mm)') = (1000 * pbp_table_retimed.('Mass (kg)') ./ (rho_water * hp_area)).*pbp_table_retimed.('SWE factor'); 
        pbp_table_retimed.('FBF SWE Accumulation (mm)') = cumsum(pbp_table_retimed.('FBF SWE (mm)'));  
        
        % calculate snow:
        
        pbp_table_retimed.('FBF Snow (mm)') = rho_water * pbp_table_retimed.('FBF SWE (mm)')./ pbp_table_retimed.('Heat Flux Density (kg/m^3)');
        pbp_table_retimed.('FBF Snow Accumulation (mm)') = cumsum(pbp_table_retimed.('FBF Snow (mm)')); 

    else 
        % create a summary table with all 0's:
        
        avi_summary_table = table(pbp_table.Time(1));
        avi_summary_table.duration = (pbp_table.Time(end)-pbp_table.Time(1)); 
        avi_summary_table.rectCx = 0;
        avi_summary_table.circleCx = 0;
        avi_summary_table.perCx = 0;
        avi_summary_table.ellipCx = 0; 
        avi_summary_table.sdi = 0;
        avi_summary_table.rho = 0;
        avi_summary_table.pbpSWE = 0;
        avi_summary_table.fbfSWE = 0;
        avi_summary_table.pbpSnow = 0;
        avi_summary_table.fbfSnow = 0;
        avi_summary_table.hotPlateArea = hp_area;
        avi_summary_table.SWEfactor = SWE_factor;
        avi_summary_table.minMass = h_mass_fbf_min;
        
        % create a blank cell for averaged time table if there are no
        % particles:
        
        pbp_table_retimed = [];
    end
    
    % add variable names:  

    avi_summary_table.Properties.VariableNames = {'Time', 'Duration', 'Rect Complexity', 'Circle Complexity', 'Perimeter Complexity', 'Ellipse Complexity', 'SDI', 'HFD Density (kg*m^-3)', 'PBP SWE (mm)', 'FBF SWE (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'Hot Plate Area', 'SWE Factor', 'Min FBF Mass'}; 
    avi_summary_table = table2timetable(avi_summary_table);

    % store as a cell so can be used outside of parforloop:

    pbp_table_cell{file_i} = pbp_table; 
    pbp_table_filtered_cell{file_i} = pbp_table_filtered;  
    avi_summary_table_cell{file_i} = avi_summary_table; 
    pbp_table_retimed_cell{file_i} = pbp_table_retimed; 
    

 end

%% now put back into one large time table 

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

% time averaged data: 

pbp_table_retimed = vertcat(pbp_table_retimed_cell{:}); 
pbp_table_retimed = sortrows(pbp_table_retimed, 'Time'); % sort by time 
pbp_table_retimed.("FBF SWE (mm)")  = fillmissing(pbp_table_retimed.("FBF SWE (mm)"),  'constant', 0);
pbp_table_retimed.("FBF Snow (mm)") = fillmissing(pbp_table_retimed.("FBF Snow (mm)"), 'constant', 0);
pbp_table_retimed.("FBF SWE Accumulation (mm)") = cumsum(pbp_table_retimed.("FBF SWE (mm)"));
pbp_table_retimed.("FBF Snow Accumulation (mm)") = cumsum(pbp_table_retimed.("FBF Snow (mm)"));

%% save processed tables

startTime = datestr(pbp_table.Time(1), 'yyyy-mm-dd_HH-MM-ss');

% unfiltered particle data table:

writetimetable(pbp_table, [output_dir, '/DEID_unfilteredParticle_', startTime, '.csv']);

% filtered particle data table:

writetimetable(pbp_table_filtered, [output_dir, '/DEID_filteredParticle_', startTime, '.csv']);

% .avi summary table:

writetimetable(avi_summary_table, [output_dir, '/DEID_aviTotals_', storm_output, '.csv']);

% time averaged data table:

writetimetable(pbp_table_retimed, [output_dir,'/DEID_TS_10min_', startTime, '.csv']);

[~, parent_dir, ~] = fileparts(pwd);
disp(['Saved Output for: ', parent_dir])

%% 

subplot(2,2,1)
histogram(pbp_table_filtered.("Rect Complexity"))
title('Rectangular Cx')
subplot(2,2,2)
histogram(pbp_table_filtered.("Circle Complexity"))
title('Circular Cx')
subplot(2,2,3)
histogram(pbp_table_filtered.("Perimeter Complexity"))
title('Perimeter Cx')
subplot(2,2,4)
histogram(pbp_table_filtered.("Ellipse Complexity"))
title('Ellipse Cx')