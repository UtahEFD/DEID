function [SWE_fbf, time_series_fbf, h_data_cells, hp_area, h_mass_fbf_min, plate_temp] = fbf_method( ...
    vid, num_frames, frame_cropped_ref, final_crop_indexes, ...
    min_thres, minimum_hydro_area, mPerPix, m2PerPix2, ...
    int_to_temp_conversion, colorbar_image_indexes, k_dLv, vid_fps, time_series)

%% "frame by frame method"; this is how we obtain SWE for each .avi file

h_data_cells = cell(num_frames,1);
sum_h_area_times_dt = nan(num_frames,1);
plate_temps = nan(num_frames,1);

for frame_ii = 1:num_frames
    % frame = frames{frame_ii};
    frame = read(vid, frame_ii);
    frame_gray = im2gray(frame); % convert frame of interest to gray scale
    frame_noCB = imcrop(frame_gray, colorbar_image_indexes);
    plate_int = double(max(frame_noCB(:)));
    plate_temp = plate_int * int_to_temp_conversion;
    plate_temps(frame_ii) = plate_temp;
    frame_cropped = imcrop(frame_gray, final_crop_indexes); % remove colorbar and kapton tape from image
    frame_filtered = frame_cropped > min_thres; % removed below min threshold, on rbg ([0, 255]) scale 
    frame_filled = imfill(frame_filtered, 'Holes'); % clean up Hydrometeors
    frame_final = bwareaopen(frame_filled, minimum_hydro_area); % any hydrometeor whose area is less than minimum_hydro_area (set to 2 pixels) is disgarded

    % now continue on to get hydrometeor properties: 

    h_geo_prop = regionprops(frame_final, 'PixelIdxList', 'MajorAxisLength', 'MinorAxisLength', 'Centroid', 'Area', 'BoundingBox', 'Perimeter'); % returns the centroid, the area , and the bounding box (left, top, width, height) of each blob
    
    % if no properties are found, go to next frame: 

    if (isempty(h_geo_prop))
        continue;
    end

    % build hydrometeor property matrices from regionprops values: 

    h_bounding_box = cat(1,h_geo_prop.BoundingBox); % concat all values to bounding box indices in pixels
    rect_widthPix = h_bounding_box(:,3); 
    rect_heightPix = h_bounding_box(:,4);
    h_centroid = round(cat(1, h_geo_prop.Centroid)); % concat all values centroid indexes
    h_perimeterPix = cat(1,h_geo_prop.Perimeter); % concat all values hydrometeor perimeters in pixels
    h_areaPix = cat(1, h_geo_prop.Area); % concat all hydrometoer areas in pixels
    h_majorPix = cat(1,h_geo_prop.MajorAxisLength); % concat all major axis in pixels
    % h_minorPix = cat(1,h_geo_prop.MinorAxisLength); % concat all minor axis in pixels
    
    % convert to length scales: 

    rect_widthM = rect_widthPix * mPerPix;
    rect_heightM = rect_heightPix * mPerPix;
    h_perimeterM = h_perimeterPix * mPerPix; 
    h_area = h_areaPix .* m2PerPix2; 
    h_majorM = h_majorPix * mPerPix;
    % h_minorM = h_minorPix * pix_to_m_conversion;

    % calculate circumscribed areas: 

    h_rectAreaM = rect_widthM .* rect_heightM; % rectangle
    % h_circleAreaM = (pi * h_majorM.^2)/4; % circumscribed circle using major axis 
    
    % difference in temperature of each centroid and the plate:
    % NOTE: centroid of an irregular blob can fall on a warm background pixel, causing dT -> 0.
    % Replaced below with per-pixel mean dT.
    % h_centroid_i = sub2ind(size(frame_cropped), h_centroid(:, 2), h_centroid(:, 1)); % find the linear index of the centriods in orginal image
    % snowflake_int = double(frame_cropped(h_centroid_i)); % intensities of centroid pixels of snow
    % plate_h_dtemp = plate_temp - (snowflake_int .* int_to_temp_conversion);

    % mean temperature difference between plate and each hydrometeor,
    % computed over all pixels in each blob (clipped to >= 0 to prevent
    % non-physical negative dT from warm centroid pixels):

    num_blobs = numel(h_geo_prop);
    plate_h_dtemp = zeros(num_blobs, 1);
    for b = 1:num_blobs
        pix_int = double(frame_cropped(h_geo_prop(b).PixelIdxList));
        pix_dT  = plate_temp - pix_int .* int_to_temp_conversion;
        plate_h_dtemp(b) = mean(max(pix_dT, 0));
    end

    % product of hydrometeor area with the mean temp difference:
    % h_area_times_dtemp = h_area .* plate_h_dtemp;  % (original centroid-based line kept for reference)
    h_area_times_dtemp = h_area .* plate_h_dtemp;
    
    % sum the product of individual area and temp. diff in each frame:
    % **this is how we obtain hydrometeor mass using fbf method**

    sum_h_area_times_dt(frame_ii)=sum(h_area_times_dtemp);         
    
    % build large matrix of Hydrometeor data:

    h_data_cells{frame_ii} = cat(2, h_centroid, plate_h_dtemp, h_perimeterM, h_area, h_rectAreaM, h_majorM); 
end

% frame by frame SWE calculation:

sum_h_area_times_dt(isnan(sum_h_area_times_dt)) =0; % turn all NaN to 0's
hp_area = numel(frame_cropped_ref) * m2PerPix2; % hotplate area - to subtract noisy areas: - mean(noisyA)        
h_mass_fbf = (k_dLv*sum_h_area_times_dt) / vid_fps; % total mass evaporates in each frame
h_mass_fbf_min = min(h_mass_fbf); % we know the plate should be empty when it is not snowing..
% h_mass_fbf = h_mass_fbf - h_mass_fbf_min; % subtract off min mass on a frame to account for any resiude
SWE_fbf = h_mass_fbf / hp_area;
time_series_fbf = time_series(1:length(SWE_fbf));
plate_temp = mean(plate_temps, 'omitnan');
end
