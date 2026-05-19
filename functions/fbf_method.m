function [SWE_fbf, time_series_fbf, h_data_cells, hp_area, h_mass_fbf_min, k_dLv_calibrate, plate_temp, back_temp] = fbf_method( ...
    vid, num_frames, frame_cropped_ref, kapton_indexes, colorbar_image_indexes, colorbar_kapton_image_indexes, ...
    min_thres, minimum_hydro_area, mPerPix, m2PerPix2, ...
    int_to_temp_conversion, k_dLv, vid_fps, time_series)

%% "frame by frame method"; this is how we obtain SWE for each .avi file

h_data_cells = cell(num_frames,1);
plate_int = nan(num_frames,1);
plate_temp = nan(num_frames,1);
back_temp = nan(num_frames,1);
q_surr = nan(num_frames,1);
k_dLv_calibrate = nan(num_frames,1);
sum_h_area_times_dt = nan(num_frames,1);

for frame_ii = 1:num_frames     
    % frame = frames{frame_ii};
    frame = read(vid, frame_ii);
    frame_gray = im2gray(frame); % convert frame of interest to gray scale
    kaptonTape = imcrop(frame_gray, kapton_indexes);% crop only the kapton tape
    plate_int(frame_ii) = max(kaptonTape(:)); % this finds max intensity of the kapton tape (plate temperature)  
    plate_temp(frame_ii) = plate_int(frame_ii)* int_to_temp_conversion; % plate temperature in degrees C 
    frame_cropped = imcrop(frame_gray, colorbar_kapton_image_indexes); % back to orginal grayscale image... now remove colorbar and kapton tape from image
    frame_filtered = frame_cropped > min_thres; % removed below min threshold, on rbg ([0, 255]) scale 
    frame_filled = imfill(frame_filtered, 'Holes'); % clean up Hydrometeors
    frame_final = bwareaopen(frame_filled, minimum_hydro_area); % any hydrometeor whose area is less than minimum_hydro_area (set to 2 pixels) is disgarded
    back_temp(frame_ii) = mean(max(double(frame_cropped)))*int_to_temp_conversion; % background temp in degrees C

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

    % calibrate k/d:
    q_surr(frame_ii) = ((5.67e-8 *(back_temp(frame_ii)*int_to_temp_conversion+273.15)^4) - (5.67e-8*0.1*(plate_temp(frame_ii)*int_to_temp_conversion+273.15)^4))/(1-0.1); 
    k_dLv_calibrate(frame_ii) = max(0, (-3.085e-5 + (8.5336e-8*min_thres*int_to_temp_conversion*plate_temp(frame_ii))...
        + (1.909e-11*(min_thres*int_to_temp_conversion*plate_temp(frame_ii) - 1393.67)^2))...
        *(q_surr(frame_ii)/min_thres));

    % q_surr(frame_ii) = ((5.67e-8 *(12+273.15)^4) - (5.67e-8*0.1*(100+273.15)^4))/(1-0.1); 
    % k_dLv_calibrate(frame_ii) = (-3.085e-5 + (8.5336e-8*min_thres*int_to_temp_conversion*100)...
    %     + (1.909e-11*((min_thres*int_to_temp_conversion*100) - 1393.67)^2))...
    %     *(q_surr(frame_ii)/min_thres); 

    % calculate circumscribed areas: 

    h_rectAreaM = rect_widthM .* rect_heightM; % rectangle
    % h_circleAreaM = (pi * h_majorM.^2)/4; % circumscribed circle using major axis 
    
    % difference in temperature of each centroid and the plate:

    h_centroid_i = sub2ind(size(frame_cropped), h_centroid(:, 2), h_centroid(:, 1)); % find the linear index of the centriods in orginal image
    snowflake_int = min(double(frame_cropped(h_centroid_i)), plate_int(frame_ii)); % intensities of centroid pixels of snow, clamped to plate intensity
    plate_h_dtemp = (plate_int(frame_ii)* int_to_temp_conversion) - (snowflake_int .* int_to_temp_conversion);
    
    % product of hydrometeor area with the temp difference:

    h_area_times_dtemp = h_area .* plate_h_dtemp;         
    
    % sum the product of individual area and temp. diff in each frame:
    % **this is how we obtain hydrometeor mass using fbf method**

    sum_h_area_times_dt(frame_ii)=sum(h_area_times_dtemp);         
    
    % build large matrix of Hydrometeor data:

    h_data_cells{frame_ii} = cat(2, h_centroid, plate_h_dtemp, h_perimeterM, h_area, h_rectAreaM, h_majorM, (ones(numel(h_geo_prop),1)*k_dLv_calibrate(frame_ii))); 
end

% frame by frame SWE calculation:

sum_h_area_times_dt(isnan(sum_h_area_times_dt)) =0; % turn all NaN to 0's
hp_area = numel(frame_cropped_ref) * m2PerPix2; % hotplate area - to subtract noisy areas: - mean(noisyA)
h_mass_fbf = (k_dLv_calibrate.*sum_h_area_times_dt) / vid_fps; % total mass evaporates in each frame
h_mass_fbf_min = min(h_mass_fbf); % we know the plate should be empty when it is not snowing..
h_mass_fbf = h_mass_fbf - h_mass_fbf_min; % subtract off min mass on a frame to account for any resiude
SWE_fbf = h_mass_fbf / hp_area;
time_series_fbf = time_series(1:length(SWE_fbf));
end
