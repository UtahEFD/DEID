# DEID (Differential Emissivity Imaging Device)
This respository includes code and documentation describing the DEID. 

The DEID is a new evaporation-based optical and thermal instrument designed to measure  mass, size, density, and type (i.e., snow, rain, and mixtures) of individual hydrometeors.

Files in this directory include:
(1) DEID_code_variables.m - this Matlab script does the following The input of this code is thermal images. The output of individual hydrometeor especially, mass, area, evaporation time, density, which can be calculated with this Matlab script. After calculation of individual parameters, SWE rate, snow rate, visibility can be estimated using the 2nd part of this code.Information about initial parameters and calibration parameter is written in the code and each line is explained.
(2) Istan_swe_rho.m this code is for instantaneous swe rate and density. every frame, 1/fps sec, swe rate can calculate. 
