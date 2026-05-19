function pbp_table = pbp_method( ...
    dT_fbf, perimeter_fbf, area_fbf, rectArea_fbf, majorAxis_fbf, max_h_obs, ...
    areaTol, k_dLv, k_dLv_calibrate, vid_fps, time_series, hp_area, mmPerM, SWE_fbf, ...
    SWEfactor_threshold, hf_rho_coeff, rho_water, rho_ice, sigma_ice)

%% "particle by particle method"

% initialize parameters:

all_h_appears_ind = []; % initial time index of each snowflake in time array
h_perimeter = {}; % hydrometeor perimeter over time
h_majorAxis = {}; % hydrometeor major axis over time
h_rectArea = {}; % circumscribed rectangle area

deltaTemp_range = []; % difference between max and min values of hydrometeor deltaTemp
deltaTemp_residue_flags = []; % flags for when deltaTemp does not change
hArea_range = []; % difference between max and min values of hydrometeor area
hArea_residue_flags = []; % flags for when area does not change 

h_max_area = {}; % hydrometeor maximum area
h_delta_time = {}; % hydrometeor time to melt and evaporate
h_delta_temp_max = {}; % hydrometeor centroid's max intensity
h_delta_temp_mean = {}; % hydrometeor centroid's mean intensity
h_mass_pbp = {}; % hydrometeor calculated mass

h_k_dLv_calibrate = []; % k_dLv calibrated

% loop over all Hydrometeors for each frame: 

for h_ii = 1:max_h_obs
    h_appear_evap_bool = diff(area_fbf(h_ii, :) > 0); % creates an array 
    % of indices that are booleans. +1 is when snowflake appears and -1
    % when it evaps. this is finding the difference in a logical array!
    h_appears_ind = find(h_appear_evap_bool > 0); % finds the index when the hydrometor appears
    
    if isempty(h_appears_ind)
        continue;
    end
    
    h_evaps_ind = find(h_appear_evap_bool < 0); % hydrometeor dissapears index

    if isempty(h_evaps_ind) % get rid of hydrometeors which do no evaporate completely
        continue; 
    else
    
    % isolate the hydrometeor properties for its lifetime: 

    h_dT_tmp = dT_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
    h_perimeter_tmp = perimeter_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
    h_area_tmp = area_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
    h_rectArea_tmp = rectArea_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
    h_majorAxis_tmp = majorAxis_fbf(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
    k_dLv_tmp = k_dLv_calibrate(h_ii, h_appears_ind(1):h_evaps_ind(end)+1);
    
    % isolate for positive values:

    h_area_tmp_bool = h_area_tmp > 0; 
    
    % finds how many 'continuous' snowflakes there are: indices of contiguous '1s':

    propstemp = regionprops(h_area_tmp_bool, 'PixelIdxList'); 

    % grab properties of each snowflake over their 'life': 

    for jj = 1:numel(propstemp)
        idx = propstemp(jj).PixelIdxList;
        all_h_appears_ind(end+1) = h_appears_ind(1) + idx(1) - 1; % global frame index of this segment's first appearance

        h_k_dLv_calibrate(end+1) = mean(k_dLv_tmp(idx));
        h_perimeter{end+1} = max(h_perimeter_tmp(idx)); % max hydrometeor perimeter over time
        h_rectArea{end+1} = max(h_rectArea_tmp(idx)); % max rectangle area over time
        h_majorAxis{end+1} = max(h_majorAxis_tmp(idx)); % major axis over time
        h_max_area{end+1} = max(h_area_tmp(idx)); % max snowflake area over time
        h_delta_time{end+1} = numel(idx); % time of snowflake's life
        h_delta_temp_max{end+1} = max(h_dT_tmp(idx)); % max temperature difference between snowflake and plate
        h_delta_temp_mean{end+1} = mean(h_dT_tmp(idx)); % mean temperature difference between snowflake and plate

        % filtering:

        deltaTemp_range(end+1) = range(h_dT_tmp(idx));
        deltaTemp_residue_flags(end+1) = (deltaTemp_range(end) == 0); % flag snowflakes whose delta temp does not change over time
        hArea_range(end+1) = range(h_area_tmp(idx));
        hArea_residue_flags(end+1) = (hArea_range(end) <= areaTol); % flag snowflakes whose area does not change over time
                        
        % multiply the area by the temperature difference for that hydrometeor:
    
        h_area_times_dT_pbp = h_area_tmp(idx) .* h_dT_tmp(idx);
    
        % now integrate (sum over snowflake's life):
    
        h_mass_pbp{end+1} = sum(h_k_dLv_calibrate(end) * h_area_times_dT_pbp); % total mass that each snoflake contains across its entire life cylce

    end

    end
end

% convert to matrixes

h_perimeter = cell2mat(h_perimeter);
h_rectArea=cell2mat(h_rectArea);
h_majorAxis = cell2mat(h_majorAxis); 

h_mass_pbp=cell2mat(h_mass_pbp); % hydrometeor mass
h_max_area=cell2mat(h_max_area); % actual area
h_delta_temp_max=cell2mat(h_delta_temp_max); % temperature diffrence between plate and water droplet using max intensity
h_delta_temp_mean=cell2mat(h_delta_temp_mean); % temperature diffrence between plate and water droplet using mean intensity 

% conversions and calculations using PBP data

h_mass_pbp = h_mass_pbp / vid_fps; 
h_dEff = ((4/pi) * h_max_area).^(1/2); % convert area to diameter of hydrometeor 
h_evap_time = cell2mat(h_delta_time) * (1 / vid_fps); % evaporation time
h_vol_sph = (3/4) * h_max_area.^(3/2); % spherical volume
h_rho_sph = h_mass_pbp ./ h_vol_sph; % density calculation: spherical assumption
% h_energy_per_time = h_mass_pbp * l_constant ./ (h_max_area .* h_evap_time); % heat flux method: energy per unit area per time
% h_height = h_evap_time .* h_delta_temp_mean; 
h_rho_hfd = (hf_rho_coeff * h_mass_pbp) ./ (h_max_area .* h_evap_time .* h_delta_temp_mean);
h_rho_hfd(~isfinite(h_rho_hfd)) = 0;
h_vol_hfd = h_mass_pbp ./ h_rho_hfd; % volume of each snowflakes using mean heat flux method density
h_vol_hfd(~isfinite(h_vol_hfd)) = 0;
h_initial_time_indices = all_h_appears_ind(1:length(h_mass_pbp));
h_initial_time = time_series(h_initial_time_indices);
eqWaterDrop_area = ((9*pi)/16)^(1/3) * (h_mass_pbp ./ rho_water).^(2/3); % surface area of equivalent water droplet (See POF Singh et al. 2023)

% Cx and SDI calculations:

complexity1 = h_rectArea ./ h_max_area; % complexity using boundingBox area (See CRST Morrison et al. 2023)
% complexity2 = h_circleArea ./ h_max_area; % complexity using circle area ((pi*majorAxis^2)/4)
% complexity3 = h_circlePerimeter ./ h_perimeter; % complexity using circle perimeter and perimeter of snowflake
% complexity4 = h_ellipseArea ./ h_max_area; % complexity using ellipse area (PCA) 

sdi = h_max_area ./ eqWaterDrop_area ; % SDI (See CRST Morrison et al. 2023)
sdi(~isfinite(sdi)) = 0;

% shear strength calculation:

shearStrength = sigma_ice .* sdi .* ((h_rho_hfd./rho_ice) .^ complexity1);
shearStrength(~isfinite(shearStrength)) = 0; 

% SWE and snow calculations:

SWE_pbp = mmPerM * h_mass_pbp ./ (rho_water * hp_area); % [mm]
SWE_factor = sum(SWE_fbf) / sum(SWE_pbp); % swe factor is used to adjust pbp SWE

if ~isfinite(SWE_factor) || SWE_factor > SWEfactor_threshold
    SWE_factor = SWEfactor_threshold; % if swe factor is abnormally high or undefined, cap it
end

SWE_fbf_particles = SWE_pbp * SWE_factor; % adjusted SWE_pbp is effectivley SWE_fbf
SWE_factor_particles = SWE_fbf_particles ./ SWE_pbp; % now calculate SWE factor for all particles
SWE_factor_particles(~isfinite(SWE_factor_particles)) = 0;

snow_pbp = rho_water * (SWE_pbp ./ h_rho_hfd); % [mm]
snow_pbp(~isfinite(snow_pbp)) = 0;
snow_fbf = rho_water * (SWE_fbf_particles ./ h_rho_hfd); % [mm]
snow_fbf(~isfinite(snow_fbf)) = 0;

% store PBP data as table for post processing:

pbp_table = table(h_initial_time', h_evap_time', ...
    h_mass_pbp', h_majorAxis', h_dEff', h_perimeter', h_max_area', h_rectArea', ...
    eqWaterDrop_area', h_rho_sph', h_rho_hfd', h_vol_hfd', h_vol_sph', ...
    h_delta_temp_max', h_delta_temp_mean', complexity1', ...
    sdi', shearStrength', SWE_pbp', SWE_fbf_particles', snow_pbp', snow_fbf', ...
    SWE_factor_particles', deltaTemp_range', hArea_range', deltaTemp_residue_flags', hArea_residue_flags', h_k_dLv_calibrate');

pbp_table.Properties.VariableNames = {'Time', 'Evap Time (s)', 'Mass (kg)', 'Major Axis (m)', 'Eff Diameter (m)', ...
    'Perimeter (m)', 'Snowflake Area (m^2)', 'Rectangle Area (m^2)', ...
    'Water Droplet Area (m^2)', 'Spherical Density (kg/m^3)', 'Heat Flux Density (kg/m^3)', ...
    'Heat Flux Volume (m^3)', 'Spherical Volume (m^3)', 'Delta Temp Max', 'Delta Temp Mean', ...
    'Complexity', 'SDI', 'Shear Strength (Pa)', 'PBP SWE (mm)', 'FBF SWE (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', ...
    'SWE factor', 'Delta Temp Range', 'Area Range', 'Delta Temp Flag', 'Area Flag', 'k/d calibrated'};

pbp_table = table2timetable(pbp_table);
pbp_table = sortrows(pbp_table, 'Time');
pbp_table.("PBP SWE Accumulation (mm)") = cumsum(pbp_table.("PBP SWE (mm)"));
pbp_table.("FBF SWE Accumulation (mm)") = cumsum(pbp_table.("FBF SWE (mm)"));
pbp_table.("PBP Snow Accumulation (mm)") = cumsum(pbp_table.("PBP Snow (mm)"));
pbp_table.("FBF Snow Accumulation (mm)") = cumsum(pbp_table.("FBF Snow (mm)"));
pbp_table.("Missing Data") = false(height(pbp_table),1);
end
