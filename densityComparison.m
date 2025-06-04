% Loads in .csv files to compare density filters. 
% 
% Ben Silberman

clear, clc
%%

workingDir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID_files/UDOT';
cd(workingDir)
%%
decTotals_fileName = 'DEID_totals_dec_2023.csv'; 
decTotals_dataFullPath = fullfile(workingDir, decTotals_fileName);
decTotals_data = readtable(decTotals_dataFullPath); 

decParticles_fileName = 'DEID_Particle_2023-12-02_19-42-22.csv'; 
decParticles_dataFullPath = fullfile(workingDir, decParticles_fileName);
decParticles_data = readtable(decParticles_dataFullPath); 

%%
plot(decParticles_data.Time, decParticles_data.DensityHFD, 'o')