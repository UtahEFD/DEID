git adfunction [DEID_storm_data_table, FULL_storm_data_table] = missingTimeFunction(particle_output_table, storm_date)
% Function to find missing time between .avi files 
% Post Processing for DEID_Processor.m

%% Parse out data for each storm

particle_output_table.Time.TimeZone = 'America/Denver'; 

if storm_date == 0105
    start_time = datetime('04-Jan-2024 11:00:00', 'TimeZone', 'America/Denver');
    end_time = datetime('05-Jan-2024 15:00:00', 'TimeZone', 'America/Denver');
    data_index = timerange(start_time, end_time, 'closed');
    DEID_storm_data_table = particle_output_table(data_index, :);

elseif storm_date == 0108
    start_time = datetime('07-Jan-2024 00:00:00', 'TimeZone', 'America/Denver');
    end_time = datetime('08-Jan-2024 06:00:00', 'TimeZone', 'America/Denver');
    data_index = timerange(start_time, end_time, 'closed');
    DEID_storm_data_table = particle_output_table(data_index, :);
end

% % jan10
% jan10_start_time = datetime('09-Jan-2024 15:00:00', 'TimeZone', 'America/Denver');
% jan10_end_time = datetime('10-Jan-2024 09:00:00', 'TimeZone', 'America/Denver');
% jan10_data_index = timerange(jan10_start_time, jan10_end_time, 'closed');
% jan10_DEID_data_table = allJan_DEID_data_table(jan10_data_index, :); 

% % feb06
% feb06_start_time = datetime('05-Feb-2024 06:00:00', 'TimeZone', 'America/Denver');
% feb06_end_time = datetime('06-Feb-2024 10:00:00', 'TimeZone', 'America/Denver');
% feb06_data_index = timerange(feb06_start_time, feb06_end_time, 'closed');
% feb06_DEID_data_table = allFeb_DEID_data_table(feb06_data_index, :);

% % feb07
% feb07_start_time = datetime('06-Feb-2024 11:00:00', 'TimeZone', 'America/Denver');
% feb07_end_time = datetime('07-Feb-2024 12:00:00', 'TimeZone', 'America/Denver');
% feb07_data_index = timerange(feb07_start_time, feb07_end_time, 'closed');
% feb07_DEID_data_table = allFeb_DEID_data_table(feb07_data_index, :);

% % feb08
% feb08_start_time = datetime('07-Feb-2024 13:00:00', 'TimeZone', 'America/Denver');
% feb08_end_time = datetime('08-Feb-2024 18:00:00', 'TimeZone', 'America/Denver');
% feb08_data_index = timerange(feb08_start_time, feb08_end_time, 'closed');
% feb08_DEID_data_table = allFeb_DEID_data_table(feb08_data_index, :);

% % feb20
% feb20_start_time = datetime('19-Feb-2024 15:00:00', 'TimeZone', 'America/Denver');
% feb20_end_time = datetime('20-Feb-2024 14:00:00', 'TimeZone', 'America/Denver');
% feb20_data_index = timerange(feb20_start_time, feb20_end_time, 'closed');
% feb20_DEID_data_table = allFeb_DEID_data_table(feb20_data_index, :);

%% Generate tables using indices of large time differences  

time_diff = diff(DEID_storm_data_table.Time); % Find difference in time column
missing_time_indices = find(time_diff > seconds(2)); % Find instances greater than 2 seconds

missing_time_indices_table = []; 

% Loop over missing time indices to build table with missing data indices: 
for i = 1:length(missing_time_indices)
    missing_time_indices_table = [missing_time_indices_table; missing_time_indices(i)];
    missing_time_indices_table = [missing_time_indices_table; missing_time_indices(i) + 1];
end

% Check that last index is included:
if ~ismember(missing_time_indices(end) +1, missing_time_indices_table);
    missing_time_indices_table = [missing_time_indices_table; missing_time_indices(end)+1];
end

% Find the average number of particles that fell per time period:
unique_times = unique(DEID_storm_data_table.Time); % Extract unique times
% times_count = arrayfun(@(x) sum(jan05_DEID_data_table.Time ==x), unique_times); % count the number of occurances for each time
% average_particles = round(mean(times_count)); % Find the average
% Consider removing this above line of code. It takes very long
average_particles = 16; 

% Create a new timetable with missing time points: 
FULL_storm_data_table = DEID_storm_data_table;

for i = 1:length(missing_time_indices)
    start_index = missing_time_indices(i);
    end_index = start_index + 1;

    start_time = DEID_storm_data_table.Time(start_index);
    end_time = DEID_storm_data_table.Time(end_index);

    % Calculate the missing number of seconds:
    missing_seconds = seconds(end_time - start_time) - 1;

    % Create the number of points needed to fill gap:
    new_entries = average_particles * missing_seconds; % how many entries are needed to fill missing time duration 
    new_times = linspace(start_time + seconds(1), end_time-seconds(1), new_entries)'; % turn this into a datetime array

    % Interpolate data to correspond with missing time:
    start_value = DEID_storm_data_table{start_index, :}; % values corresponding to last row of .avi file
    end_value = DEID_storm_data_table{end_index, :}; % values corresponding to first row of next .avi file 

    new_values = arrayfun(@(k) linspace(start_value(k), end_value(k), length(new_times))', 1:size(DEID_storm_data_table, 2), 'UniformOutput', false);
    new_values = cell2mat(new_values);

    % Place new data into filled time table
    new_row = array2timetable(new_values, 'RowTimes', new_times, 'VariableNames', DEID_storm_data_table.Properties.VariableNames);
    FULL_storm_data_table = [FULL_storm_data_table; new_row];
end

% Recalculate accumulated SWE and Snow per storm
FULL_storm_data_table.SWE_Acc_mm = cumsum(FULL_storm_data_table.SWE_mm);
FULL_storm_data_table.Snow_Acc_mm = cumsum(FULL_storm_data_table.Snow_mm);
FULL_storm_data_table.Snow_Acc_in = FULL_storm_data_table.Snow_Acc_mm * (1/25.4); 

% Sort new table 
FULL_storm_data_table = sortrows(FULL_storm_data_table); 