# DEID (Differential Emissivity Imaging Device)
This respository includes the DEID code base and associated documentation. 
The DEID is a new evaporation-based optical and thermal instrument designed to measure  mass, size, density, and type (i.e., snow, rain, and mixtures) of individual hydrometeors.
This code base is used to extract hydrometeor properties from the DEID thermal camera video. Analysis of output data is done independently.

## **DEID_Processor.m**
The primary Matlab script that derives hydrometeor and bulk properties. Takes in DEID video files in AVI format and outputs both a time series of bulk hydrometeor properties as well as each individual particle recorded by DEID. 

#### Time Series Output:
- Timestamp
- SWE (total)
- Snow (total)
- Density (mean)
- Complexity (mean)
- SDI (mean)

#### Particle Table Output:
- Timestamp
- SWE 
- Snow
- Density
- Mass
- Volume
- Complexity
- SDI

## **sortPositions_v2.m**
Helper function for DEID_Processor script to sort positions of hydrometeors on hotplate surface.

## **terminalVelocity.m**
Helper function for DEID_Processor script to computer terminal velocity of each hydrometeor. 
