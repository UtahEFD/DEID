# DEID (Differential Emissivity Imaging Device)
This respository includes the DEID code base and associated documentation. 
The DEID is a new evaporation-based optical and thermal instrument designed to measure  mass, size, density, and type (i.e., snow, rain, and mixtures) of individual hydrometeors.
This code base is used to extract hydrometeor properties from the DEID thermal camera video. Analysis of output data is done independently.

## **DEID_Processor.m**
The primary Matlab script that derives hydrometeor and bulk properties. Takes in DEID video files in AVI format and outputs both a time series of bulk hydrometeor properties as well as each individual particle recorded by DEID. 

#### Time Series Output:
The time series takes an input 'time_interval' to average / sum particle values. Default set to 300 seconds. 
- Timestamp
- Terminal Velocity (mean)
- Complexity (mean)
- SDI (mean)
- Density (mean)
- SWE (total)
- Snow (total)

#### Particle Table Output:
Measures all paramaters that fall on DEID hotplate across a given .avi file. 
- Timestamp
- Terminal Velocity
- Complexity
- SDI
- Mass
- Volume
- Density
- SWE
- Snow
- SWE 
- Snow

#### Diagnostic Table Output:
Used to check outputs being calculated in code. This can be adjusted to check any desired output. Currently checks:
- SWE factor
- Number of particles

## **DEID_Processor_v2.m**
Modified version of DEID_Processor.m with modification detailed below:

#### Start/End Video Time Table:
- Creates table 'start_end_time_table' which includes all video files within directory, their start times, end times, and length of video
- Specify start and end time of a storm to index directory

#### Fills gaps between .avi files
- Calls missingTimeFunction.m which does the following:
  - finds the differences between every particle's timestamps
  - locates all the times this is greater than 60 seconds
  - multiples the time difference by average number of particles per second to create a range to average on
  - takes the average of this number of rows 

#### Time Table Output:
- Timestamp
- Terminal Velocity (mean)
- Complexity (mean)
- SDI (mean) 
- Mass (total) 
- Volume (total)
- Density (mean)
- Effective Diameter (mean)
- Surface Area (mean)
- Void Space (mean)
- Temperature Difference (plate and snowflake)
- SWE (total)
- Snow (total)

#### Particle Table Output:
- Timestamp
- Terminal Velocity
- Complexity
- SDI
- Mass
- Volume
- Density
- Effective Diameter
- Surface Area
- Void Space
- Temperature Difference (plate and snowflake)
- SWE
- Snow 

## **sortPositions_v2.m**
Helper function for DEID_Processor script to sort positions of hydrometeors on hotplate surface.

## **terminalVelocity.m**
Helper function for DEID_Processor script to computer terminal velocity of each hydrometeor. 
