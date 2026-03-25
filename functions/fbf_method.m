function [SWE_fbf, time_series_fbf, h_data_cells, hp_area, h_mass_fbf_min] = fbf_method( ...
    vid, num_frames, frame_cropped_ref, colorbar_image_indexes, colorbar_kapton_image_indexes, ...
    min_thres, minimum_hydro_area, mPerPix, m2PerPix2, ...
    int_to_temp_conversion, k_dLv, vid_fps, time_series)

%% "frame by frame method"; this is how we obtain SWE for each .avi file

h_data_cells = cell(num_frames,1);
plate_int = nan(num_frames,1);
sum_h_area_times_dt = nan(num_frames,1);

for frame_ii = 1:num_frames     
    % frame = frames{frame_ii};
    frame = read(vid, frame_ii);
    frame_gray = im2gray(frame); % convert frame of interest to gray scale
    frame_gray_cropped_wKapton = imcrop(frame_gray, colorbar_image_indexes);% crop out colorbar
    plate_int(frame_ii) = max(max(double(frame_gray_cropped_wKapton))); % this assumes max temperature in image is the plate temperature with Kapton tape 
    frame_cropped = imcrop(frame_gray, colorbar_kapton_image_indexes); % back to orginal grayscale image... now remove colorbar and kapton tape from image
    frame_filtered = frame_cropped > min_thres; % removed below min threshold, on rbg ([0, 255]) scale 
    frame_filled = imfill(frame_filtered, 'Holes'); % clean up Hydrometeors
    frame_final = bwareaopen(frame_filled, minimum_hydro_area); % any hydrometeor whose area is less than minimum_hydro_area (set to 2 pixels) is disgarded
    
    % remove centroids that appear more than 1000 times:
    % 
    % props = regionprops(frame_final, 'Area', 'Centroid','PixelIdxList');
    % 
    % if ~isempty(props)
    % 
    %     frame_centroids = cat(1, props.Centroid); % collect centroids of particles on frame 
    %     frame_area = cat(1, props.Area); % collect areas of particles on frame 
    % 
    % if isempty(noiseCentroids)
    % 
    %     isNoise = false(size(frame_centroids,1),1); % basically do nothing
    % 
    % else
    % 
    %     isNoise = ismember(frame_centroids, noiseCentroids, 'rows'); % logical array identifying which centroids on the frame are noisy ones 
    % 
    % end
    % 
    %     noisyA(frame_ii) = sum(frame_area(isNoise, :)); % store area to subtract from hpArea later 
    % 
    %     % black out noisy hydrometeors:
    % 
    %     for k = find(isNoise)'
    %         frame_final(props(k).PixelIdxList) = 0;
    %     end
    % 
    % end

    % now continue on to get hydrometeor properties: 

    h_geo_prop = regionprops(frame_final, 'PixelIdxList', 'MajorAxisLength', 'MinorAxisLength', 'Centroid', 'Area', 'BoundingBox', 'Perimeter'); % returns the centroid, the area , and the bounding box (left, top, width, height) of each blob
    
    % if no properties are found, go to next frame: 

    if (isempty(h_geo_prop))
        continue;
    end

    % PCA-BASED CIRCUMSCRIBED ELLIPSE AREA PER HYDROMETEOR
    % 
    % h_PCAellipseAreaM = zeros(length(h_geo_prop),1);
    % 
    % for ii = 1:length(h_geo_prop)
    % 
    %     % extract pixel coordinates of hydrometeor:
    % 
    %     pixList = h_geo_prop(ii).PixelIdxList;
    %     [r, c] = ind2sub(size(frame_final), pixList);
    %     pts = [c, r];  % nx2 array of pixel coordinates for each hydrometeor 
    % 
    %     % center the points:
    %     C    = mean(pts,1);   % centroid of hydrometeor in pixel coordinates 
    %     pts0 = pts - C;       % centered coordinates of pixels 
    % 
    %     % covariance + eigen decomposition:
    %     Sigma    = cov(double(pts0));
    %     [V, D]   = eig(Sigma);
    % 
    %     % sort eigenvectors to ensure major & minor axes:
    %     [~, idx] = sort(diag(D), 'descend');
    %     V = V(:, idx);
    %     D = diag(sort(diag(D), 'descend'));
    % 
    %     % rotate points into PCA basis:
    %     ptsRot = pts0 * V;
    % 
    %     % first: maximum absolute extent in PCA axes (inscribed semi-axes):
    %     a0_pix = max(abs(ptsRot(:,1)));   % initial semi-major axis
    %     b0_pix = max(abs(ptsRot(:,2)));   % initial semi-minor axis
    % 
    %     % scale ellipse so it CIRCUMSCRIBES all points:
    %     normVals = (ptsRot(:,1)/a0_pix).^2 + (ptsRot(:,2)/b0_pix).^2;
    %     s        = sqrt(max(normVals));   % >= 1
    % 
    %     % final circumscribing semi-axes (keep same variable names):
    %     a_pix = s * a0_pix + 0.5;
    %     b_pix = s * b0_pix + 0.5;
    % 
    %     % area of ellipse (pixel units) using circumscribing semi-axes:
    %     ellipse_area_pix = pi * a_pix * b_pix;
    % 
    %     % convert to m^2:
    %     h_PCAellipseAreaM(ii) = ellipse_area_pix * pix_to_m2_conversion;
    % 
    %     % checking complexity:
    %     areaPix = numel(pixList); 
    %     Cx_pix = ellipse_area_pix / areaPix;
    %     if Cx_pix < 1
    %         fprintf('Frame %d, hydro %d: Cx_pix = %.3f\n', frame_ii, ii, Cx_pix);
    %     end
    % 
    % end

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

    h_centroid_i = sub2ind(size(frame_cropped), h_centroid(:, 2), h_centroid(:, 1)); % find the linear index of the centriods in orginal image
    snowflake_int = double(frame_cropped(h_centroid_i)); % intensities of centroid pixels of snow
    plate_h_dtemp = (plate_int(frame_ii)* int_to_temp_conversion) - (snowflake_int .* int_to_temp_conversion); 
    
    % product of hydrometeor area with the temp difference:

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
h_mass_fbf = h_mass_fbf - h_mass_fbf_min; % subtract off min mass on a frame to account for any resiude
SWE_fbf = h_mass_fbf / hp_area;
time_series_fbf = time_series(1:length(SWE_fbf));
end
