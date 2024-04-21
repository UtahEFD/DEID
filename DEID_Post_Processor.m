%% DEID Post Processing Code
% Outputs data derived from particle by particle csv file
% AUTHOR : Dhiraj Singh, Benjamin Silberman, Alex Blackmer
% LAST UPDATED: 04/10/2024
clear, clc, close all
close all
clear

workingDir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID/Atwater/JAN/Jan_09_10_storm/processed_output_data';
cd(workingDir)

% Defines global variables
time_step = minutes(1);        % Resampling period in minutes
A_hot = 0.0101;       % Area of hot plate?
mm_to_inches = 1/25.4;

%% Handles FBF SWE data
fbf_table_raw = readtable('01_09_24_1349_50_FBF.csv');
fbf_table_raw.Properties.VariableNames{'Var1'} = 'Timestamp';
fbf_table_raw.Timestamp = datetime(fbf_table_raw.Timestamp);
fbf_table_raw = table2timetable(fbf_table_raw);
% Resamples time series at desired interval
fbf_table = retime(fbf_table_raw, 'regular', 'sum', 'TimeStep', time_step);

fbf_table.SWE_accum_mm = cumsum(fbf_table.SWE_fbf);
fbf_table.SWE_accum_in = fbf_table.SWE_accum_mm * mm_to_inches;

%% Handles PBP Data
% Reads in particle by particle table data
pbp_table_raw=readtable('01_09_24_1349_50_PBP.csv');
% Sorts table by timestamp
timestamp = datetime(pbp_table_raw.Hydrometeor_initial_times);
pbp_table_raw = sortrows(pbp_table_raw, 'Hydrometeor_initial_times');
mass = pbp_table_raw.Hydrometeor_Mass_pbp;    
diameter = pbp_table_raw.Hydrometeor_Diameter;    
max_area = pbp_table_raw.Hydrometeor_Max_Area;    
max_circ_area = pbp_table_raw.Hydrometeor_Max_Circumscribed_Area; 
time_to_evap = pbp_table_raw.Time_2_evap;        
delta_T = pbp_table_raw.Hydrometeor_Delta_T;    
delta_T1 = pbp_table_raw.Hydrometeor_Delta_T1; 
density_sph = pbp_table_raw.Sphere_Density;    
% Reads in heat flux density directly 
hfd_1 = pbp_table_raw.HeatFlux_Density1;
hfd_2 = pbp_table_raw.HeatFlux_Density2;                % Unused

% Computes terminal velocity
% Constants
rho_air = 0.9;
g = 9.8;
pi = 3.14;
eta = 1.81*10^-5;
% Derived values
Area_ratio = max_circ_area./max_area;
nu = eta/rho_air;
x1 = 8*mass.*g*rho_air;
x2 = pi*eta^2;
x3 = Area_ratio.^0.25;
X = x1.*x3./x2;
p1 = (1+0.1519.*X.^0.5).^0.5;
Re = 8.5*(p1-1).^2;
v = Re.*eta.*(pi./max_circ_area).^0.5;
v_t = v./(2.*rho_air);
vvol = mass./1000;
A_m = 1.2*vvol.^(2/3);
cpx1 = max_circ_area./max_area;
sdi = max_area./A_m;
mu = 1.5*10^-5;
v_st = (1/(18*mu))* hfd_1.*diameter.^2;

% Stores derived data in table for easy access
pbp_table_raw = table(timestamp, mass, density_sph, hfd_1, v_t, max_area, cpx1, sdi, diameter, max_circ_area, time_to_evap);
pbp_table_raw = sortrows(pbp_table_raw, 'timestamp');
% Filters data to find where 0 < mass < .005
[g1,g2] = find(pbp_table_raw.mass > 0 & pbp_table_raw.mass < 0.005);
pbp_table_raw = pbp_table_raw(g1,:);

% Appends volume data to table 
pbp_table_raw.VV1 = pbp_table_raw.mass./(pbp_table_raw.density_sph);
pbp_table_raw.VV2 = pbp_table_raw.mass./(pbp_table_raw.hfd_1);
pbp_table_raw = table2timetable(pbp_table_raw);

%% Computes averaged and summed PBP data
% Averages PBP data 
avg_cols = {'density_sph', 'hfd_1', 'v_t','cpx1', 'sdi', 'diameter'};
avg_table = retime(pbp_table_raw(:, avg_cols), 'regular', 'mean', 'TimeStep', time_step);
% Sums PBP data
sum_cols = {'mass', 'max_area', 'max_circ_area', 'VV1', 'VV2'};
sum_table = retime(pbp_table_raw(:, sum_cols), 'regular', 'sum', 'TimeStep', time_step);

pbp_table = horzcat(avg_table, sum_table);

%% Gets average SWE from averaged PBP data
pbp_table.SWE = pbp_table.mass ./ (A_hot);       % mm/min rate
pbp_table.SWE_accum_mm = cumsum(pbp_table.SWE);  % in of swe
pbp_table.SWE_accum_in = pbp_table.SWE_accum_mm * mm_to_inches;  % in of swe
% Finds difference factor between FBF swe1 and PBP swe2 and adjusts swe2
factor = fbf_table.SWE_accum_in(end) / pbp_table.SWE_accum_in(end);
pbp_table.SWE_F = pbp_table.SWE * factor;                                   
pbp_table.SWE_F_accum_mm = cumsum(pbp_table.SWE_F);            % mm/hr rate
pbp_table.SWE_F_accum_in = pbp_table.SWE_F_accum_mm * mm_to_inches; % in of swe

%% Gets Snow from PBP data
factor2 = mean(pbp_table.VV1) ./ mean(pbp_table.VV2);
V3 = factor2 * pbp_table.VV2;
% ddensity1 = avg_table.total_mass ./ avg_table.rho_h;      % Not used       
ddensity2 = pbp_table.mass ./ pbp_table.VV2;   % kg/m^3
% den_sphere=mean(avg_data(:,3))                            % not sure what this is
pbp_table.snow = 1000 * pbp_table.SWE_F ./ ddensity2;                      % mm/min rate
pbp_table.snow_acc_mm = cumsum(pbp_table.snow);        % mm of snow 
pbp_table.snow_acc_in = pbp_table.snow_acc_mm * mm_to_inches;  % in of snow

% %% Outputs to file
% Data_FBF_AVG = table(hd, swe2_accum_rate_f,swe2_accumulation_f,...
%     swe2_accumulation_inch_f, snow_accum_rate_ave, ...
%     snow_accumulation_ave, snow_accumulation_ave_inch, ddensity2);
% Data_FBF_AVG.Properties.VariableNames = {'Time_5_min_avg', 'SWE_Rate_mm/hr','SWE_Accum_mm',...
%     'SWE_Accum_in',  'Snow_Rate_mm/hr', 'Snow_Accum_mm', 'Snow_Accum_in', 'Density'};   
% writetable(Data_FBF_AVG, 'test2_FBF.csv');
% 
