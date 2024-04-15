%% DEID Post Processing Code
% Outputs data derived from particle by particle csv file
% AUTHOR : Dhiraj Singh, Benjamin Silberman, Alex Blackmer
% LAST UPDATED: 04/10/2024
clear, clc, close all
close all
clear all

workingDir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID/Atwater/JAN/Jan_09_10_storm/processed_output_data/';
cd(workingDir)

% Defines global variables
A_hot=0.0056;       % Area of hot plate?
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
swe1_rate=  swe1_table.SWE_fbf * frame_rate * one_hour_in_seconds;  % SWE rate in mm/hr
swe1_time = 0:del:del*(length(swe1_table.SWE_fbf)-1);% Number of seconds
nod1 = avg_min * 60;                                 % Number of data indices.
num_periods = floor(length(swe1_table.SWE_fbf)/nod1);% Number of avg periods
swe1_time_series = [swe1_time' swe1_rate];
swe1_period_series = zeros(num_periods,2);
% Fills swe1 period series with averaged values from time series
for j = 1:num_periods  
    swe1_period_series(j,:) = mean(swe1_time_series(nod1*(j-1)+1:nod1*(j-1)+nod1,:));
end

time_hour = nod1 / (frame_rate * one_hour_in_seconds); 
swe1_accu = swe1_period_series(:,2) * (time_hour);
swe1_accumulation = cumsum(swe1_accu);
swe1_accumulation_inch = swe1_accumulation / 25.4;
% Manipulates time series to give date and time rather than just seconds
t_time = (t_c + swe1_period_series(:,1));
tt_time = 27475200+864000 + t_time;% variable

%% Handles PBP Data
% Reads in particle by particle table data
pbp_table=readtable('01_09_24_1349_50_PBP.csv');
h_initial_time_index = pbp_table.Hydrometeor_initial_time_indexes;
h_initial_time_index = pbp_table.Hydrometeor_initial_time_indexes;
m = pbp_table.Hydrometeor_Mass_pbp;    
D = pbp_table.Hydrometeor_Diameter;    
A = pbp_table.Hydrometeor_Max_Area;    
AA = pbp_table.Hydrometeor_Max_Circumscribed_Area; 
t_evp = pbp_table.Time_2_evap;        
d_sph = pbp_table.Sphere_Density;     % Why is d_sph duplicated?
TT = pbp_table.Hydrometeor_Delta_T;    
TT1 = pbp_table.Hydrometeor_Delta_T1; 

d_heat = 6.4418e+04*m./(A.*t_evp.*TT);
d_sph = 6.4418e+04*m./(A.*t_evp.*TT1);  % Why is d_sph duplicated?
raw_den_heat = mean(d_heat);

% Computes terminal velocity
% Constants
rho_air = 0.9;
g = 9.8;
pi = 3.14;
eta = 1.81*10^-5;
% Derived values
Area_ratio = AA./A;
nu = eta/rho_air;
x1 = 8*m.*g*rho_air;
x2 = pi*eta^2;
x3 = Area_ratio.^0.25;
X = x1.*x3./x2;
p1 = (1+0.1519.*X.^0.5).^0.5;
Re = 8.5*(p1-1).^2;
v = Re.*eta.*(pi./AA).^0.5;
v_t = v./(2.*rho_air);
vvol = m./1000;
A_m = 1.2*vvol.^(2/3);
cpx1 = AA./A;
sdi = A./A_m;
mu = 1.5*10^-5;
v_st = (1/(18*mu))* d_heat.*D.^2;
% Stores derived data in array
data_pbp = [h_initial_time_index m d_sph d_heat v_t A cpx1 sdi D AA t_evp];
% Filters data to find where 0 < mass < .005
[g1,g2] = find(data_pbp(:,2)>0 & data_pbp(:,2)<0.005);
data_pbp = data_pbp(g1,:);
% Appends velocity data to data array
VV1 = data_pbp(:,2)./(data_pbp(:,3));
VV2 = data_pbp(:,2)./(data_pbp(:,4));
data_pbp = [data_pbp VV1 VV2];

%% Averages PBP data
for i = 1:one_hour_in_seconds
    ini = (data_pbp(:,1)>=(one_hour_in_minutes*(i-1)+0));
    fin = (data_pbp(:,1)<=(one_hour_in_minutes*(i-1)+one_hour_in_minutes));
    avg_time(i) = nanmean(data_pbp(ini & fin,1)); % mean time
    total_mass(i) = sum(data_pbp(ini & fin,2)); % total mass in given time
    rho_s(i) = nanmean(data_pbp(ini & fin,3)); % mean sph. density in given time
    rho_h(i) = nanmean(data_pbp(ini & fin,4)); % mean heat density in given time
    avg_velocity(i) = nanmean(data_pbp(ini & fin,5));
    total_area(i) = sum(data_pbp(ini & fin,6)); % mean density in 1 min
    avg_cpx(i) = nanmean(data_pbp(ini & fin,7));    
    avg_sdi(i) = nanmean(data_pbp(ini & fin,8));  
    avg_dia(i) = nanmean(data_pbp(ini & fin,9)); 
    total_cir_area(i) = sum(data_pbp(ini & fin,10));
    total_vol1(i) = sum(data_pbp(ini & fin,11));
    total_vol2(i) = sum(data_pbp(ini & fin,12));
end
% Stores final data and removes missing values
avg_data = [avg_time' total_mass' rho_s' rho_h'  avg_velocity' total_area' avg_cpx' avg_sdi' avg_dia' total_cir_area' total_vol1' total_vol2' ];
avg_data = rmmissing(avg_data);

%% Gets average SWE from averaged PBP data
r_time = (t_c + avg_data(:,1));
rr_time = 27475200+864000+r_time;
hd = datenum(seconds(rr_time));

swe2 = seconds_per_min*avg_data(:,2)./(A_hot);          % mm/min rate
swe2_accum_rate = swe2.*(1/seconds_per_min);            % mm/hr rate
swe2_accumulation_inch = cumsum(swe2_accum_rate) * mm_to_inches;  % in of swe
% Finds difference factor between swe1 and swe2 and adjusts swe2
factor = swe1_accumulation_inch(end)/swe2_accumulation_inch(end);
swe2_f = swe2*factor;                                   
swe2_accum_rate_f = swe2_f.*(1/seconds_per_min);        % mm/min rate
swe2_accumulation_f = cumsum(swe2_accum_rate_f);        % mm/hr rate
swe2_accumulation_inch_f = swe2_accumulation_f * mm_to_inches; % in of swe

%% Gets Snow from PBP data
factor2 = mean(avg_data(:,11))./mean(avg_data(:,12));
V3 = factor2*avg_data(:,12);
ddensity1 = avg_data(:,2)./avg_data(:,4);               % kg/m^3 ?
ddensity2 = avg_data(:,2)./avg_data(:,12);
% den_sphere=mean(B(:,3)) % not sure what this is
snow_ave = 1000*swe2_f./ddensity1;                      % mm/min rate
snow_accum_rate_ave = snow_ave.*(1/seconds_per_min);    % mm/hr rate
snow_accumulation_ave = cumsum(snow_accum_rate_ave);    % mm of snow 
snow_accumulation_ave_inch = snow_accumulation_ave * mm_to_inches;  % in of snow

%% Outputs to file
Data_FBF_AVG = table(hd, swe2_accum_rate_f,swe2_accumulation_f,...
    swe2_accumulation_inch_f, snow_accum_rate_ave, ...
    snow_accumulation_ave, snow_accumulation_ave_inch, ddensity1);
Data_FBF_AVG.Properties.VariableNames = {'Time_5_min_avg', 'SWE_Rate_mm/hr','SWE_Accum_mm',...
    'SWE_Accum_in',  'Snow_Rate_mm/hr', 'Snow_Accum_mm', 'Snow_Accum_in', 'Density'};   
writetable(Data_FBF_AVG, 'test_FBF.csv');

