function deid = get_deid_params()

deid.colorbar_image_indexes = [1 1 383 288];

deid.final_crop_indexes = [30 30 340 240];

% calibration constants
deid.k_dLv = 0.003;

% Singh et al. (2024)
deid.hf_rho_coeff = 1.01e05; % [K*s*m^-1]

end