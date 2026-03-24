function [dT_fbf, perimeter_fbf, area_fbf, rectArea_fbf, majorAxis_fbf, h_data_sorted, max_h_obs] = ...
    sort_h_data_cells(h_data_cells, num_frames, sort_threshold, filename)

%% call sortPositions_v2.m

num_cols = [];
for ii = 1:num_frames
    if ~isempty(h_data_cells{ii})
        num_cols = size(h_data_cells{ii}, 2);
        break
    end
end

if isempty(num_cols)
    warning('Video %s contains no hydrometeors. Skipping.', filename);
    dT_fbf = [];
    perimeter_fbf = [];
    area_fbf = [];
    rectArea_fbf = [];
    majorAxis_fbf = [];
    h_data_sorted = {};
    max_h_obs = [];
    return
end

for ii = 1:num_frames
    if isempty(h_data_cells{ii})
        h_data_cells{ii} = zeros(0, num_cols);
    end
end

h_data_sorted = cell(size(h_data_cells));
h_data_sorted{1} = h_data_cells{1};

for frame_jj = 2:num_frames
    h_data_sorted{frame_jj} = h_data_cells{frame_jj};
    h_data_sorted{frame_jj} = sortPositions_v2(h_data_sorted{frame_jj-1}, h_data_sorted{frame_jj}, sort_threshold);
end

% return frame with max number of hydrometeors:

max_h_obs = max(cellfun(@(x) size(x, 1), h_data_sorted, 'UniformOutput', 1));

% now pad the data with zeros so the frames all have the same number:

h_data_sorted = cellfun(@(x) cat(1, x, zeros(max_h_obs - size(x, 1), width(h_data_sorted{1}))), ...
    h_data_sorted, 'UniformOutput', 0);

%% isolating the variables and put them into a matrix to work with:

% for reference: [h_centroid(1), h_centroid(2), plate_h_dtemp,... 
% h_perimeterM, h_area, h_rectAreaM, h_majorM]

dT_fbf = cellfun(@(x) x(:, 3), h_data_sorted, 'UniformOutput', 0);
perimeter_fbf = cellfun(@(x) x(:, 4), h_data_sorted, 'UniformOutput', 0);
area_fbf = cellfun(@(x) x(:, 5), h_data_sorted, 'UniformOutput', 0);
rectArea_fbf = cellfun(@(x) x(:, 6), h_data_sorted, 'UniformOutput', 0);
majorAxis_fbf = cellfun(@(x) x(:, 7), h_data_sorted, 'UniformOutput', 0);

dT_fbf = cat(2,dT_fbf{:});
perimeter_fbf = cat(2, perimeter_fbf{:});
area_fbf = cat(2,area_fbf{:});
rectArea_fbf = cat(2,rectArea_fbf{:});
majorAxis_fbf = cat(2,majorAxis_fbf{:});
end
