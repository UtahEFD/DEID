function [file_names, vid_date, storm_output] = get_sorted_videos(working_dir)
%GET_SORTED_VIDEOS Find and sort .avi files by timestamp
%
% Inputs:
%   working_dir - path to directory containing video files
%
% Outputs:
%   file_names  - cell array of sorted .avi filenames
%   vid_date    - datetime array of corresponding file timestamps
%   storm_output - formatted date string (DD-MM-YYYY) from first video

    % Validate input
    if ~isfolder(working_dir)
        error('Provided working_dir is not a valid folder.');
    end

    % Get directory contents
    directory = dir(working_dir);

    % Preallocate (max possible size)
    file_names = cell(1, length(directory));
    vid_date   = NaT(1, length(directory)); % datetime array

    count = 0;

    % Loop through directory items
    for file_i = 1:length(directory)
        name = directory(file_i).name;

        if ~directory(file_i).isdir && endsWith(name, '.avi', 'IgnoreCase', true)
            count = count + 1;

            file_names{count} = name;

            % Use full path for safety
            full_path = fullfile(working_dir, name);
            vid_info = dir(full_path);

            vid_date(count) = datetime(vid_info.date);
        end
    end

    % Trim unused entries
    file_names = file_names(1:count);
    vid_date   = vid_date(1:count);

    % Handle empty case
    if isempty(file_names)
        storm_output = '';
        warning('No .avi files found in the directory.');
        return;
    end

    % Format first date
    storm_output = datestr(vid_date(1), 'dd-mm-yyyy');

    % Sort by date
    [vid_date, sort_idx] = sort(vid_date);
    file_names = file_names(sort_idx);

end