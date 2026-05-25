function phys = get_physical_constants(sensor)

% UNIT CONVERSIONS

if strcmp(sensor, 'snowflake3')
    phys.mPerPix = 3.1750e-04;
elseif strcmp(sensor, 'snowflake4')
    phys.mPerPix = 0.01/40;
else
    error('get_physical_constants: unknown sensor "%s". Expected "snowflake3" or "snowflake4".', sensor);
end

phys.m2PerPix2 = phys.mPerPix^2;
phys.mmPerM = 1e3;

% temperature/intensity relationship
phys.max_temp = 145;
phys.max_int = 255;
phys.int_to_temp_conversion = phys.max_temp / phys.max_int;

% PHYSICAL CONSTANTS

phys.rho_water = 1000;     % [kg/m^3]
phys.mu = 1.5e-5;         % [kg/(m*s)]
phys.rho_ice = 917;          % [kg/m^3]
phys.sigma_ice = 124;        % Pa 

end
