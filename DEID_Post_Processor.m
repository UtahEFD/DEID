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
A_hot = 0.0101;       % Area of hot plate?
avg_min = 5;        % Average period in minutes
frame_rate = 15;
t_c=  2340;         % Starting time?
one_hour_in_seconds = 3600;    % One hour in seconds
one_hour_in_minutes = 60;
seconds_per_min = 60;
mm_to_inches = 1/25.4;

%% Handles FBF SWE data
% swe1 is derived from frame by frame method
swe1_table = readtable('01_09_24_1349_50_FBF.csv');
del = 1/frame_rate;    % Data frequency in seconds. Inverse of frame rate
swe1_rate =  swe1_table.SWE_fbf * frame_rate * one_hour_in_seconds;  % SWE rate in mm/hr
swe1_time = 0:del:del*(length(swe1_table.SWE_fbf)-1);% Number of seconds
nod1 = avg_min * 60;                                 % Average period in seconds.
num_periods = floor(length(swe1_table.SWE_fbf)/nod1);% Number of avg periods
swe1_time_series = [swe1_time' swe1_rate];
swe1_period_series = zeros(num_periods,2);
% Fills swe1 period series with averaged values from time series
for j = 1:num_periods  
    swe1_period_series(j,:) = mean(swe1_time_series(nod1*(j-1)+1:nod1*(j-1)+nod1,:));
end

time_hour = nod1 / (frame_rate * one_hour_in_seconds); 
swe1_accu = swe1_period_series(:,2) * (time_hour);      % mm/hr
swe1_accumulation = cumsum(swe1_accu);
swe1_accumulation_inch = swe1_accumulation / 25.4;
% Manipulates time series to give date and time rather than just seconds
t_time = (t_c + swe1_period_series(:,1));
tt_time = 27475200+864000 + t_time;% variable

%% Testing new time averaging
% Reads in particle by particle table data
pbp_table2=readtable('01_09_24_1349_50_PBP.csv');


%% Handles PBP Data
% Reads in particle by particle table data
pbp_table=readtable('01_09_24_1349_50_PBP.csv');
h_initial_time_index = pbp_table.Hydrometeor_initial_time_indexes;
mass = pbp_table.Hydrometeor_Mass_pbp;    
diameter = pbp_table.Hydrometeor_Diameter;    
max_area = pbp_table.Hydrometeor_Max_Area;    
max_circ_area = pbp_table.Hydrometeor_Max_Circumscribed_Area; 
time_to_evap = pbp_table.Time_2_evap;        
delta_T = pbp_table.Hydrometeor_Delta_T;    
delta_T1 = pbp_table.Hydrometeor_Delta_T1; 
density_sph = pbp_table.Sphere_Density;    
% Reads in heat flux density directly 
hfd_1 = pbp_table.HeatFlux_Density1;
hfd_2 = pbp_table.HeatFlux_Density2;                % Unused

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
pbp_table = table(h_initial_time_index, mass, density_sph, hfd_1, v_t, max_area, cpx1, sdi, diameter, max_circ_area, time_to_evap);
% Filters data to find where 0 < mass < .005
[g1,g2] = find(pbp_table.mass > 0 & pbp_table.mass < 0.005);
pbp_table = pbp_table(g1,:);

% Appends volume data to table 
pbp_table.VV1 = pbp_table.mass./(pbp_table.density_sph);
pbp_table.VV2 = pbp_table.mass./(pbp_table.hfd_1);

%% Averages PBP data
% Define column names
column_names = {'avg_time', 'total_mass', 'rho_s', 'rho_h', 'avg_vel', ...
    'total_area', 'avg_cpx', 'avg_sdi', 'avg_dia', 'total_cir_area', ...
    'total_vol1', 'total_vol2' }; 
avg_table = array2table(nan(one_hour_in_seconds, numel(column_names)), 'VariableNames', column_names);

parfor i = 1:one_hour_in_seconds
    % Aggregates and averages particles found in time period
    ini = (pbp_table.h_initial_time_index >= (one_hour_in_minutes * (i-1) + 0));
    fin = (pbp_table.h_initial_time_index <= (one_hour_in_minutes * (i-1) + one_hour_in_minutes));
    logical_index = ini & fin;
    fd = pbp_table(logical_index,:);    % Filtered table
    avg_table(i).avg_time = mean(fd.h_initial_time_index, 'omitnan'); % mean time
    avg_table(i).total_mass = sum(fd.mass, 'omitnan'); % total mass in given time
    avg_table(i).rho_s = mean(fd.density_sph, 'omitnan'); % mean sph. density in given time
    avg_table(i).rho_h = mean(fd.hfd_1, 'omitnan'); % mean heat density in given time
    avg_table(i).avg_vel = mean(fd.v_t, 'omitnan');
    avg_table(i).total_area = sum(fd.max_area, 'omitnan'); % mean density in 1 min
    avg_table(i).avg_cpx = mean(fd.cpx1, 'omitnan');    
    avg_table(i).avg_sdi = mean(fd.sdi, 'omitnan');  
    avg_table(i).avg_dia = mean(fd.diameter, 'omitnan'); 
    avg_table(i).total_cir_area = sum(fd.max_circ_area, 'omitnan');
    avg_table(i).total_vol1 = sum(fd.VV1, 'omitnan');
    avg_table(i).total_vol2 = sum(fd.VV2, 'omitnan');
end

% Removes missing values where no snow was measured
avg_table = rmmissing(avg_table);

%% Gets average SWE from averaged PBP data
r_time = (t_c + avg_table.avg_time);
rr_time = 27475200+864000+r_time;
hd = datenum(seconds(rr_time));

swe2 = seconds_per_min * avg_table.total_mass ./ (A_hot);       % mm/min rate
swe2_accum_rate = swe2 .* (1/seconds_per_min);                  % mm/hr rate
swe2_accumulation_inch = cumsum(swe2_accum_rate) * mm_to_inches;  % in of swe
% Finds difference factor between FBF swe1 and PBP swe2 and adjusts swe2
factor = swe1_accumulation_inch(end) / swe2_accumulation_inch(end);
swe2_f = swe2 * factor;                                   
swe2_accum_rate_f = swe2_f .* (1 / seconds_per_min);        % mm/min rate
swe2_accumulation_f = cumsum(swe2_accum_rate_f);            % mm/hr rate
swe2_accumulation_inch_f = swe2_accumulation_f * mm_to_inches; % in of swe

%% Gets Snow from PBP data
factor2 = mean(avg_table.total_vol1) ./ mean(avg_table.total_vol2);
V3 = factor2 * avg_table.total_vol2;
% ddensity1 = avg_table.total_mass ./ avg_table.rho_h;      % Not used       
ddensity2 = avg_table.total_mass ./ avg_table.total_vol2;   % kg/m^3
% den_sphere=mean(avg_data(:,3))                            % not sure what this is
snow_ave = 1000 * swe2_f ./ ddensity2;                      % mm/min rate
snow_accum_rate_ave = snow_ave .* (1/seconds_per_min);      % mm/hr rate
snow_accumulation_ave = cumsum(snow_accum_rate_ave);        % mm of snow 
snow_accumulation_ave_inch = snow_accumulation_ave * mm_to_inches;  % in of snow

%% Outputs to file
Data_FBF_AVG = table(hd, swe2_accum_rate_f,swe2_accumulation_f,...
    swe2_accumulation_inch_f, snow_accum_rate_ave, ...
    snow_accumulation_ave, snow_accumulation_ave_inch, ddensity2);
Data_FBF_AVG.Properties.VariableNames = {'Time_5_min_avg', 'SWE_Rate_mm/hr','SWE_Accum_mm',...
    'SWE_Accum_in',  'Snow_Rate_mm/hr', 'Snow_Accum_mm', 'Snow_Accum_in', 'Density'};   
writetable(Data_FBF_AVG, 'test_FBF.csv');

