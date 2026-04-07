function [pbp_table, avi_summary_table, pbp_table_filtered] = append_gap_row_and_summary( ...
    pbp_table, hp_area, rho_water, rho_ice, sigma_ice, h_mass_fbf_min, evapTime_min, evapTime_max)

% Fill a 1-minute gap row with summary data and build AVI summary

gapFillWindow = minutes(1);

if height(pbp_table) > 0
    pbp_table_filtered = pbp_table( ...
        pbp_table.('Area Flag') ~= 1 & ...
        pbp_table.('Delta Temp Flag') ~= 1 & ...
        pbp_table.("Evap Time (s)") > evapTime_min & ...
        pbp_table.("Evap Time (s)") < evapTime_max, :);
else
    pbp_table_filtered = timetable();
end

% if height(pbp_table_filtered) > 0
%     final_time = pbp_table_filtered.Time(end);
%     prev_time = final_time - gapFillWindow;
%     prev_data = pbp_table_filtered( ...
%         pbp_table_filtered.Time >= prev_time & ...
%         pbp_table_filtered.Time <= final_time, :);

if height(pbp_table) > 0
    final_time = pbp_table.Time(end);
    prev_time = final_time - gapFillWindow;
    prev_data = pbp_table( ...
        pbp_table.Time >= prev_time & ...
        pbp_table.Time <= final_time, :);

    timeRow = final_time + seconds(30);

    evapTimeRow = mean(prev_data.("Evap Time (s)"));
    massRow = sum(prev_data.("Mass (kg)"));
    majAxisRow = mean(prev_data.("Major Axis (m)"));
    effDiaRow = mean(prev_data.("Eff Diameter (m)"));
    perRow = mean(prev_data.("Perimeter (m)"));
    areaRow = mean(prev_data.("Snowflake Area (m^2)"));
    rectAreaRow = mean(prev_data.("Rectangle Area (m^2)"));
    waterDropletAreaRow = mean(prev_data.("Water Droplet Area (m^2)"));
    volumeHFDrow = sum(prev_data.("Heat Flux Volume (m^3)"));
    volumeSPHrow = sum(prev_data.("Spherical Volume (m^3)"));
    deltaTempMaxRow = mean(prev_data.("Delta Temp Max"));
    deltaTempMeanRow = mean(prev_data.("Delta Temp Mean"));
    cx1Row = mean(prev_data.("Complexity"));
    sdiRow = mean(prev_data.SDI);
    SWEfactorRow = mean(prev_data.("SWE factor"));

    tempRangeRow = NaN;
    aRangeRow = NaN;
    tempFlagRow = NaN;
    aFlagRow = NaN;

    densityHFDrow = massRow / volumeHFDrow;
    densitySPHrow = massRow / volumeSPHrow;
    shearStrengthRow = sigma_ice * sdiRow * ((densityHFDrow / rho_ice) ^ cx1Row);

    PBPsweRow = 1000 * massRow / (rho_water * hp_area);
    FBFsweRow = PBPsweRow * SWEfactorRow;
    PBPsnowRow = rho_water * (PBPsweRow ./ densityHFDrow);
    FBFsnowRow = rho_water * (FBFsweRow ./ densityHFDrow);

    new_row = table(timeRow, ...
        evapTimeRow, massRow, majAxisRow, effDiaRow, perRow, areaRow, rectAreaRow, waterDropletAreaRow, ...
        densitySPHrow, densityHFDrow, volumeHFDrow, volumeSPHrow, deltaTempMaxRow, ...
        deltaTempMeanRow, cx1Row, sdiRow, shearStrengthRow, PBPsweRow, FBFsweRow, ...
        sum(pbp_table.("PBP SWE (mm)")) + PBPsweRow, ...
        sum(pbp_table.("FBF SWE (mm)")) + FBFsweRow, ...
        PBPsnowRow, FBFsnowRow, ...
        sum(pbp_table.("PBP Snow (mm)")) + PBPsnowRow, ...
        sum(pbp_table.("FBF Snow (mm)")) + FBFsnowRow, SWEfactorRow, ...
        tempRangeRow, aRangeRow, tempFlagRow, aFlagRow, ...
        'VariableNames', {'Time', 'Evap Time (s)', 'Mass (kg)', 'Major Axis (m)', 'Eff Diameter (m)', ...
        'Perimeter (m)', 'Snowflake Area (m^2)', 'Rectangle Area (m^2)', ...
        'Water Droplet Area (m^2)', 'Spherical Density (kg/m^3)', 'Heat Flux Density (kg/m^3)', ...
        'Heat Flux Volume (m^3)', 'Spherical Volume (m^3)', 'Delta Temp Max', 'Delta Temp Mean', ...
        'Complexity', 'SDI', 'Shear Strength (Pa)', 'PBP SWE (mm)', 'FBF SWE (mm)', ...
        'PBP SWE Accumulation (mm)', 'FBF SWE Accumulation (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', ...
        'PBP Snow Accumulation (mm)', 'FBF Snow Accumulation (mm)', 'SWE factor', ...
        'Delta Temp Range', 'Area Range', 'Delta Temp Flag', 'Area Flag'});

    new_row.('Missing Data') = true;
    new_row = table2timetable(new_row);

    pbp_table = [pbp_table; new_row];
    pbp_table_filtered = [pbp_table_filtered; new_row];

    pbp_table_filtered.("PBP SWE Accumulation (mm)") = cumsum(pbp_table_filtered.("PBP SWE (mm)"));
    pbp_table_filtered.("FBF SWE Accumulation (mm)") = cumsum(pbp_table_filtered.("FBF SWE (mm)"));
    pbp_table_filtered.("PBP Snow Accumulation (mm)") = cumsum(pbp_table_filtered.("PBP Snow (mm)"));
    pbp_table_filtered.("FBF Snow Accumulation (mm)") = cumsum(pbp_table_filtered.("FBF Snow (mm)"));

    avi_summary_table = build_avi_summary_table( ...
        pbp_table_filtered, hp_area, SWEfactorRow, h_mass_fbf_min, rho_ice, sigma_ice);

else
    avi_summary_table = table(NaT);
    avi_summary_table.duration = seconds(NaN);
    avi_summary_table.Cx = 0;
    avi_summary_table.sdi = 0;
    avi_summary_table.rho = 0;
    avi_summary_table.shearStrength = 0;
    avi_summary_table.pbpSWE = 0;
    avi_summary_table.fbfSWE = 0;
    avi_summary_table.pbpSnow = 0;
    avi_summary_table.fbfSnow = 0;
    avi_summary_table.hotPlateArea = hp_area;
    avi_summary_table.SWEfactor = 0;
    avi_summary_table.minMass = h_mass_fbf_min;
    avi_summary_table.Properties.VariableNames = {'Time', 'Duration', 'Complexity', 'SDI', 'Density (kg*m^-3)', ...
        'Shear Strength (Pa)', 'PBP SWE (mm)', 'FBF SWE (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', ...
        'Hot Plate Area', 'SWE Factor', 'Min FBF Mass'};
    avi_summary_table = table2timetable(avi_summary_table);
end

end