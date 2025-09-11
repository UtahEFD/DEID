% Loads in .csv files to compare density filters. 
% 
% Ben Silberman

clear, clc
%% load files 

% DEID directory: 
deidDataPath = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/Atwater/test';
pbp_file = 'DEID_Particle_2023-03-05_17-59-09.csv'; 
summary_file = 'DEID_totals_mar07_2023.csv'; 

% read in files as tables: 
deidTimeTable_pbp = readtable(fullfile(deidDataPath,pbp_file),'VariableNamingRule','preserve');
deidTimeTable_sum = readtable(fullfile(deidDataPath,summary_file),'VariableNamingRule','preserve');



%%
figure()
plot(deidTimeTable_pbp.Time, deidTimeTable_pbp.("Density Sphere"), 'o')
title('Density (kg m^{-3}) Time Series')
xlabel('Time')
ylabel('Density (kg m^{-3})')
%%
figure()
plot(deidTimeTable_pbp.Time, deidTimeTable_pbp.("Evap Time"), 'o')
title('Evap Time (s) Time Series')
xlabel('Time')
ylabel('Evap Time (s)')

figure()
plot(deidTimeTable_pbp.Time, deidTimeTable_pbp.("Temp Diff"), 'o')
title('Temp Diff (s) Time Series')
xlabel('Time')
ylabel('Temp Diff')
