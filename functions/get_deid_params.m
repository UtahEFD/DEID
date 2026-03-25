function deid = get_deid_params()

deid.colorbar_image_indexes = [1 1 384 288];

deid.crop_index = 43;

deid.colorbar_kapton_image_indexes = [ ...
    1, ...
    (deid.colorbar_image_indexes(2) + deid.crop_index), ...
    383, ...
    (deid.colorbar_image_indexes(4) - deid.crop_index) ...
];

% calibration constants
deid.k_dLv = 0.002;

% Singh et al. (2024)
deid.hf_rho_coeff = 1.01e05; % [K*s*m^-1]

end