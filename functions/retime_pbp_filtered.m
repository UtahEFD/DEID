function pbp_table_retimed = retime_pbp_filtered(pbp_table_filtered, time_step, rho_water, hp_area)

% RETIME_PBP_FILTERED  Retime a filtered particle timetable to a regular grid
%
% pbp_table_retimed = retime_pbp_filtered(pbp_table_filtered, time_step, rho_water, hp_area)
%
% Inputs:
%   pbp_table_filtered  - timetable (or table with variable 'Time') of per-particle rows
%   time_step           - duration, e.g. minutes(10)
%   rho_water           - numeric (e.g. 1000)
%   hp_area             - hot plate area (same units used in your pipeline)
%
% Output:
%   pbp_table_retimed   - timetable with regular bins at 'time_step' and computed:
%       - mean columns for avg_cols
%       - sum columns for sum_cols
%       - Density (kg/m^3), FBF/PBP SWE (mm), FBF/PBP Snow (mm)
%       - Accumulation columns (cumulative sums) computed by treating NaN
%         as zero *for contributions only* (so cumsum stays flat across gaps).
%
% Notes:
%  - preserves NaN for averaged columns in empty bins (so gaps remain gaps)
%  - does not permanently overwrite the original re-timed sum/mean columns
%  - input may be a table with a 'Time' column or a timetable
%
% Example:
%  pbp_table_retimed = retime_pbp_filtered(pbp_table_filtered, minutes(10), 1000, 0.005);

% -----------------------
% Input checks / quick exits
% -----------------------

if nargin < 4
    error('All four inputs required: pbp_table_filtered, time_step, rho_water, hp_area');
end

% If empty input, return empty timetable

if isempty(pbp_table_filtered)
    pbp_table_retimed = timetable();
    return;
end

% If input is table, convert to timetable (expecting a 'Time' variable)

if ~istimetable(pbp_table_filtered)
    if ismember('Time', pbp_table_filtered.Properties.VariableNames)
        pbp_table_filtered = table2timetable(pbp_table_filtered, 'RowTimes', pbp_table_filtered.Time);
    else
        error('Input must be a timetable or a table with a ''Time'' variable.');
    end
end

% Ensure time_step is a duration:

if ~isduration(time_step)
    error('time_step must be a duration, e.g. minutes(10)');
end

% -----------------------
% Define columns
% -----------------------

avg_cols = {'Complexity','SDI','Eff Diameter (m)', 'Major Axis (m)', 'Snowflake Area (m^2)','Evap Time (s)','SWE factor'};
sum_cols = {'Mass (kg)','Heat Flux Volume (m^3)'};

% Add missing columns as NaN so retime doesn't fail (keeps schema consistent)
allCols = [avg_cols, sum_cols];
for ii = 1:numel(allCols)
    cn = allCols{ii};
    if ~ismember(cn, pbp_table_filtered.Properties.VariableNames)
        pbp_table_filtered.(cn) = NaN(height(pbp_table_filtered),1);
    end
end

% -----------------------
% Build time grid automatically from min->max (regular bins)
% -----------------------

allTimes = pbp_table_filtered.Properties.RowTimes;
tmin = min(allTimes);
tmax = max(allTimes);

% Construct regular grid covering [tmin, tmax]
% Let MATLAB create bins automatically by using 'regular' with TimeStep (simple method)
% (If you need custom anchoring, change this section)

avg_table = retime(pbp_table_filtered(:, avg_cols), 'regular', 'mean', 'TimeStep', time_step);
sum_table = retime(pbp_table_filtered(:, sum_cols), 'regular', 'sum', 'TimeStep', time_step);

pbp_table_retimed = horzcat(avg_table, sum_table);

% -----------------------
% Derived variables (no temp columns)
% -----------------------

% Mass and heat volume (may be NaN in empty bins)
mass = pbp_table_retimed.('Mass (kg)');
heatVol = pbp_table_retimed.('Heat Flux Volume (m^3)');

% Density: NaN unless both mass and heatVol are present and heatVol != 0
density = NaN(size(mass));
validMask = ~isnan(mass) & ~isnan(heatVol) & (heatVol ~= 0);
density(validMask) = mass(validMask) ./ heatVol(validMask);
pbp_table_retimed.('Density (kg/m^3)') = density;

% Ensure SWE factor present
if ~ismember('SWE factor', pbp_table_retimed.Properties.VariableNames)
    pbp_table_retimed.('SWE factor') = NaN(height(pbp_table_retimed),1);
end
swe_factor = pbp_table_retimed.('SWE factor');

% PBP / FBF SWE per bin (NaN if missing inputs)
PBP_SWE = NaN(size(mass));
FBF_SWE = NaN(size(mass));
validMass = ~isnan(mass) & ~isnan(swe_factor) & (hp_area ~= 0);
if any(validMass)
    PBP_SWE(validMass) = (1000 .* mass(validMass) ./ (rho_water * hp_area));
    FBF_SWE(validMass) = PBP_SWE(validMass) .* swe_factor(validMass);
end
pbp_table_retimed.('PBP SWE (mm)') = PBP_SWE;
pbp_table_retimed.('FBF SWE (mm)') = FBF_SWE;

% -----------------------
% Accumulations: treat NaN as 0 only for contribution (do not overwrite original columns)
% -----------------------
PBP_SWE_contrib = PBP_SWE; PBP_SWE_contrib(isnan(PBP_SWE_contrib)) = 0;
FBF_SWE_contrib = FBF_SWE; FBF_SWE_contrib(isnan(FBF_SWE_contrib)) = 0;

pbp_table_retimed.('PBP SWE Accumulation (mm)') = cumsum(PBP_SWE_contrib);
pbp_table_retimed.('FBF SWE Accumulation (mm)') = cumsum(FBF_SWE_contrib);

% Snow per bin (requires density)
FBF_Snow = NaN(size(FBF_SWE));
PBP_Snow = NaN(size(PBP_SWE));
validFBF = ~isnan(FBF_SWE) & ~isnan(density) & (density ~= 0);
validPBP = ~isnan(PBP_SWE) & ~isnan(density) & (density ~= 0);

FBF_Snow(validFBF) = rho_water .* FBF_SWE(validFBF) ./ density(validFBF);
PBP_Snow(validPBP) = rho_water .* PBP_SWE(validPBP) ./ density(validPBP);

pbp_table_retimed.('FBF Snow (mm)') = FBF_Snow;
pbp_table_retimed.('PBP Snow (mm)') = PBP_Snow;

% Snow accumulations (NaN -> 0 for contribution)
FBF_Snow_contrib = FBF_Snow; FBF_Snow_contrib(isnan(FBF_Snow_contrib)) = 0;
PBP_Snow_contrib = PBP_Snow; PBP_Snow_contrib(isnan(PBP_Snow_contrib)) = 0;

pbp_table_retimed.('FBF Snow Accumulation (mm)') = cumsum(FBF_Snow_contrib);
pbp_table_retimed.('PBP Snow Accumulation (mm)') = cumsum(PBP_Snow_contrib);

end