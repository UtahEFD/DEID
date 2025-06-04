% Loads in .csv files to compare density filters. 
% 
% Ben Silberman

clear, clc
%% load files 

% DEID directory: 
deidDataPath = '/Volumes/benji/modelData/deidData';
% DEID files:

% pbp:
deiDataFile_pbp_storm1 = 'DEID_Particle_2023-03-05_17-59-09.csv';
% deidDataFile_pbp_storm2 = 
% deidDataFile_pbp_storm3 =  
deidDataFile_pbp_storm4 = 'DEID_Particle_2023-03-19_20-02-32.csv';
deidDataFile_pbp_storm5 = 'DEID_Particle_2023-03-24_01-04-18.csv';
deidDataFile_pbp_storm6 = 'DEID_Particle_2023-12-02_19-42-22.csv'; 
deidDataFile_pbp_storm7 = 'DEID_Particle_2023-12-08_00-34-34.csv'; 
deidDataFile_pbp_storm8 = 'DEID_Particle_2024-01-04_09-33-36.csv';
deidDataFile_pbp_storm9 = 'DEID_Particle_2024-01-06_17-07-42.csv';
deidDataFile_pbp_storm10 = 'DEID_Particle_2024-01-09_13-51-00.csv';
deidDataFile_pbp_storm11 = 'DEID_Particle_2024-02-05_07-06-35.csv';
deidDataFile_pbp_storm12 = 'DEID_Particle_2024-02-06_17-04-01.csv';
deidDataFile_pbp_storm13 = 'DEID_Particle_2024-02-07_12-39-54.csv';
deidDataFile_pbp_storm14 = 'DEID_Particle_2024-02-19_14-59-28.csv';
% deidDataFile_pbp_storm15 = '';
deidDataFile_pbp_storm16 = 'DEID_Particle_2024-12-29_15-41-22.csv';
deidDataFile_pbp_storm17 = 'DEID_Particle_2025-01-01_15-45-14.csv';
deidDataFile_pbp_storm18 = 'DEID_Particle_2025-02-13_08-44-19.csv';
deidDataFile_pbp_storm19 = 'DEID_Particle_2025-03-31_15-05-01.csv';

% summary:
deidDataFile_total_storm1 = 'DEID_totals_mar0723.csv';
% deidDataFile_total_storm2
% deidDataFile_total_storm3
deidDataFile_total_storm4 = 'DEID_totals_mar2023.csv';
deidDataFile_total_storm5 = 'DEID_totals_mar2423.csv';
deidDataFile_total_storm6 = 'DEID_totals_dec0423.csv';
deidDataFile_total_storm7 = 'DEID_totals_dec0923.csv';
deidDataFile_total_storm8 = 'DEID_totals_jan0524.csv';
deidDataFile_total_storm9 = 'DEID_totals_jan0824.csv';
deidDataFile_total_storm10 = 'DEID_totals_jan1024.csv';
deidDataFile_total_storm11 = 'DEID_totals_feb0624.csv';
deidDataFile_total_storm12 = 'DEID_totals_feb0724.csv';
deidDataFile_total_storm13 = 'DEID_totals_feb0824.csv';
deidDataFile_total_storm14 = 'DEID_totals_feb2024.csv';
% deidDataFile_total_storm15 
deidDataFile_total_storm16 = 'DEID_totals_dec2924.csv';
deidDataFile_total_storm17 = 'DEID_totals_jan0125.csv';
deidDataFile_total_storm18 = 'DEID_totals_feb1325.csv';
deidDataFile_total_storm19 = 'DEID_totals_mar3125.csv';

%% read files and store as a cell array 

% place all deid pbp files into an array:   
deidDataFiles_pbp_allStorms = {deiDataFile_pbp_storm1, deidDataFile_pbp_storm4, ...
    deidDataFile_pbp_storm5, deidDataFile_pbp_storm6, deidDataFile_pbp_storm7, ...
    deidDataFile_pbp_storm8, deidDataFile_pbp_storm9, deidDataFile_pbp_storm10, ...
    deidDataFile_pbp_storm11, deidDataFile_pbp_storm12, deidDataFile_pbp_storm13, ...
    deidDataFile_pbp_storm14, deidDataFile_pbp_storm16, deidDataFile_pbp_storm17, ...
    deidDataFile_pbp_storm18, deidDataFile_pbp_storm19};

% place all deid summary files into an array: 
deidDataFiles_total_allStorms = {deidDataFile_total_storm1, deidDataFile_total_storm4, ...
    deidDataFile_total_storm5, deidDataFile_total_storm6, deidDataFile_total_storm7, ...
    deidDataFile_total_storm8, deidDataFile_total_storm9, deidDataFile_total_storm10, ...
    deidDataFile_total_storm11, deidDataFile_total_storm12, deidDataFile_total_storm13, ...
    deidDataFile_total_storm14, deidDataFile_total_storm16, deidDataFile_total_storm17, ...
    deidDataFile_total_storm18, deidDataFile_total_storm19};


deidTimeTable_pbp = cell(1,length(deidDataFiles_pbp_allStorms));
deidTimeTable_sum = cell(1,length(deidDataFiles_total_allStorms));
for i = 1:length(deidDataFiles_pbp_allStorms)
    deidTimeTable_pbp{i} = readtable(fullfile(deidDataPath,deidDataFiles_pbp_allStorms{i}),'VariableNamingRule','preserve');
    deidTimeTable_sum{i} = readtable(fullfile(deidDataPath,deidDataFiles_total_allStorms{i}),'VariableNamingRule','preserve');
end

%%

plot(deidTimeTable_pbp{1}.Time, deidTimeTable_pbp{1}.("Density HFD"), 'o')
hold on
plot(deidTimeTable_sum{1}.Time, deidTimeTable_sum{1}.HFDDensity_kg_m__3_, 'o')
hold off
legend('pbp', 'averages')