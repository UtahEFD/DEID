function [full_storm_data_table] = missingTimeFunction(particle_output_table)
% Function to find missing time between .avi files 
% Post Processing for DEID_Processor.m 

%% Generate a table using indices of time differences greater than 60 seconds:  
timeDiff = diff(particle_output_table.Time); % Find difference in time column
missingTimeIndices = find(timeDiff > seconds(60)); % Find instances greater than 2 seconds

% this can be adjusted, however seems to be accurate from testing a few
% storms:
% gather the average number of particles per .avi file:
averageParticles = 16; 

%% Create a new timetable with missing time points: 
for i = 1:length(missingTimeIndices)

    startIndex = missingTimeIndices(i);
    endIndex = startIndex + 1;
    startTime = particle_output_table.Time(startIndex);
    endTime = particle_output_table.Time(endIndex);

    % Calculate the missing number of seconds:
    missingSeconds = seconds(endTime - startTime) - 1;

    % Find the number of rows to average on:
    particlesToAverage = averageParticles * missingSeconds; % how many previous entries to average 
    
    % Take the average of these rows:
    newParticleValues = mean(particle_output_table((startIndex-round(particlesToAverage)):startIndex, :));
    newParticleValues.Time = NaT;

    % compile exising data with new data point 
    if i == 1
        particleOutputTableTemp = particle_output_table(i:startIndex, :);
        full_storm_data_table = [particleOutputTableTemp;newParticleValues];
    else
        particleOutputTableTemp = particle_output_table(missingTimeIndices(i-1)+1:startIndex, :);
        full_storm_data_table = [full_storm_data_table;particleOutputTableTemp;newParticleValues];   
    end
end

% place remaining data into new table:
full_storm_data_table = [full_storm_data_table; particle_output_table(missingTimeIndices(end)+1:end, :)]; 

% convert to timetable:
full_storm_data_table = table2timetable(full_storm_data_table);