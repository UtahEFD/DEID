function [pbp_table, pbp_table_filtered, avi_summary_table] = process_one_video( ...
    filename, working_dir, ...
    colorbar_image_indexes, colorbar_kapton_image_indexes, ...
    min_thres, minimum_hydro_area, sort_threshold, areaTol, ...
    SWEfactor_threshold, evapTime_min, evapTime_max, ...
    mPerPix, m2PerPix2, int_to_temp_conversion, ...
    k_dLv, hf_rho_coeff, rho_water, mmPerM)

% PROCESS_ONE_VIDEO Process one AVI file end-to-end.

full_filename = fullfile(working_dir, filename);
vid = VideoReader(full_filename);

% define a stable cropped image for sizes (and optional intensity
% fallback):

tmp = readFrame(vid);
tmpg = im2gray(tmp);
frame_cropped_ref = imcrop(tmpg, colorbar_kapton_image_indexes);

% restart reader so you can read from the beginning again:

vid = VideoReader(full_filename);

% get metaData for video processing:

vid_dir = dir(full_filename);
vid_fps = vid.FrameRate;
num_frames = vid.NumFrames;

% construct time series:

vid_end_time = datetime([vid_dir.date]);
time_series = datetime(vid_end_time - (0:num_frames) * seconds(1/vid_fps), 'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
time_series = flip(time_series);

% CALL fbf_method.m and run frame-by-frame method: this outputs SWE_fbf

[SWE_fbf, ~, h_data_cells, hp_area, h_mass_fbf_min] = fbf_method( ...
    vid, num_frames, frame_cropped_ref, colorbar_image_indexes, colorbar_kapton_image_indexes, ...
    min_thres, minimum_hydro_area, mPerPix, m2PerPix2, ...
    int_to_temp_conversion, k_dLv, vid_fps, time_series);

% CALL sort_h_data_cells.m organizes data and calls sortPositions_v2.m

[dT_fbf, perimeter_fbf, area_fbf, rectArea_fbf, majorAxis_fbf, ~, max_h_obs] = ...
    sort_h_data_cells(h_data_cells, num_frames, sort_threshold, filename);

if isempty(max_h_obs)
    pbp_table = timetable();
    pbp_table_filtered = timetable();
    avi_summary_table = timetable();
    return;
end

% CALL pbp_method.m and run particle by particle method:

pbp_table = pbp_method( ...
    dT_fbf, perimeter_fbf, area_fbf, rectArea_fbf, majorAxis_fbf, max_h_obs, ...
    areaTol, k_dLv, vid_fps, time_series, hp_area, mmPerM, SWE_fbf, ...
    SWEfactor_threshold, hf_rho_coeff, rho_water);

% FILTER data: call append_gap_row_and_summary for filtered pbp data
% append_gap_row_and_summary will call build_avi_summary_table.m

[pbp_table, avi_summary_table, pbp_table_filtered] = append_gap_row_and_summary( ...
    pbp_table, hp_area, rho_water, h_mass_fbf_min, evapTime_min, evapTime_max);
end
