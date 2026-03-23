% loads filteredParticle data and aviTotals to preform time averaging
clear, clc

% directories and paths:

working_dir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/processedData/202502';
output_dir = working_dir;
filteredFile = dir(fullfile(working_dir, 'DEID_filteredParticle*.csv'));
totalsFile = dir(fullfile(working_dir, 'DEID_aviTotals*.csv'));

% load files:

avi_summary_table = readtimetable(fullfile(working_dir, totalsFile(1).name), VariableNamingRule="preserve");
pbp_table_filtered = readtimetable(fullfile(working_dir, filteredFile(1).name), VariableNamingRule="preserve");

%%  set variables:

time_step = seconds(600);
rho_water = 1000; % [kg/m^3]
hp_area = avi_summary_table.("Hot Plate Area")(1); % m^2

% call function to retime

pbp_table_retimed = retime_pbp_filtered(pbp_table_filtered, time_step, rho_water, hp_area);

% save TS output: 

saveTime = datestr(pbp_table_retimed.Time(1), 'yyyy-mm-dd_HH-MM-ss');
writetimetable(pbp_table_retimed, [output_dir,'/DEID_TS_10min_', saveTime, '.csv']);

%% quick plots:

plot(pbp_table_retimed.Time, pbp_table_retimed.("FBF SWE Accumulation (mm)"))
