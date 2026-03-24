function avi_summary_table = build_avi_summary_table(pbp_table_filtered, hp_area, SWE_factor, h_mass_fbf_min)

%BUILD_AVI_SUMMARY_TABLE Create summary row for one AVI file.

avi_summary_table = table(pbp_table_filtered.Time(1));
avi_summary_table.duration = (pbp_table_filtered.Time(end)-pbp_table_filtered.Time(1));
avi_summary_table.rectCx = mean(pbp_table_filtered.('Complexity'));
avi_summary_table.sdi = mean(pbp_table_filtered.SDI);
avi_summary_table.rho = sum(pbp_table_filtered.("Mass (kg)")) / sum(pbp_table_filtered.("Heat Flux Volume (m^3)"));
avi_summary_table.pbpSWE = pbp_table_filtered.("PBP SWE Accumulation (mm)")(end);
avi_summary_table.fbfSWE = pbp_table_filtered.("FBF SWE Accumulation (mm)")(end);
avi_summary_table.pbpSnow = pbp_table_filtered.("PBP Snow Accumulation (mm)")(end);
avi_summary_table.fbfSnow = pbp_table_filtered.("FBF Snow Accumulation (mm)")(end);
avi_summary_table.hotPlateArea = hp_area;
avi_summary_table.SWEfactor = SWE_factor;
avi_summary_table.minMass = h_mass_fbf_min;
avi_summary_table.Properties.VariableNames = {'Time', 'Duration', 'Complexity', 'SDI', 'Density (kg*m^-3)', 'PBP SWE (mm)', 'FBF SWE (mm)', 'PBP Snow (mm)', 'FBF Snow (mm)', 'Hot Plate Area', 'SWE Factor', 'Min FBF Mass'};
avi_summary_table = table2timetable(avi_summary_table);
end
