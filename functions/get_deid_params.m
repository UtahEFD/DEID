function deid = get_deid_params()

deid.kapton_indexes = [190 1 50 10];
deid.colorbar_image_indexes = [1 1 384 288];

deid.crop_index = deid.kapton_indexes(4)+5;

deid.colorbar_kapton_image_indexes = [ ...
    1+ deid.crop_index, ...
    (1 + deid.crop_index), ...
    deid.colorbar_image_indexes(3) - deid.crop_index, ...
    (deid.colorbar_image_indexes(4) - deid.crop_index) ...
];

% calibration constants
deid.k_dLv = 0.003;

% Singh et al. (2024)
deid.hf_rho_coeff = 1.01e05; % [K*s*m^-1]

end