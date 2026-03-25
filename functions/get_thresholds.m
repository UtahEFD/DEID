function thresh = get_thresholds()

thresh.min_thres = 70;
thresh.minimum_hydro_area = 2;
thresh.sort_threshold = 20;
thresh.minimum_drop_life = 0;

thresh.areaTol = 0;
thresh.SWEfactor_threshold = 1.85;

thresh.evapTime_min = 1/15;
thresh.evapTime_max = 30;

thresh.noiseThresh = 999;

end