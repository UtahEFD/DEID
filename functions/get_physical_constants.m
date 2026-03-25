function phys = get_physical_constants()

% UNIT CONVERSIONS

% phys.mPerPix = .01/40; % (2025 - present)
phys.mPerPix = 3.1750e-04; % (2023 - 2024)

phys.m2PerPix2 = phys.mPerPix^2;
phys.mmPerM = 1e3;

% temperature/intensity relationship
phys.max_temp = 145;
phys.max_int = 255;
phys.int_to_temp_conversion = phys.max_temp / phys.max_int;

% PHYSICAL CONSTANTS

phys.rho_water = 1000;     % [kg/m^3]
phys.mu = 1.5e-5;         % [kg/(m*s)]

end
