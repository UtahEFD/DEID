clearvars;
close all
clc
%% DEID AVI File Processing Code
% Use to process input avi files to output txt files
% This code is based on Dhiraj_DEID_Code.m
% Authors: Dhiraj Singh, Alex Blackmer, Ben Silberman
% Date: 4/3/2024

%% Sets path and defines directories: 
workingDir = '/uufs/chpc.utah.edu/common/home/snowflake3/DEID/';
inputDir = 'Atwater/Jan_09_10_storm/';
outputDir = 'processed_output_data/';
% Sets path for .m to run in:
cd([workingDir, inputDir]) 

% Check if the output directory exists
if ~exist(outputDir, 'dir')
    % If the directory doesn't exist, create it
    mkdir(outputDir);
    disp(['Directory "', outputDir, '" created successfully.']);
else
    % If the directory already exists, display a message
    disp(['Directory "', outputDir, '" already exists.']);
end

%% Identifies all .avi files in inputDir 
% Get a list of all files and folders in this folder
filesAndFolders = dir(".");
% Initialize an empty cell array to store the file names
fileNames = {};
% Loop through each item in the input directory
for i = 1:length(filesAndFolders)
    % Check if the item is a file (not a folder) and if it ends with .avi
    if ~filesAndFolders(i).isdir && endsWith(filesAndFolders(i).name, '.avi', 'IgnoreCase', true)
        % Get the name of the file and append it to the list
        fileNames{end+1} = filesAndFolders(i).name;
    end
end

%% Processes each .avi file present in input directory
for file_i = 1:length(fileNames)
    disp(['Processing File: ', fileNames{file_i}])
    vid=VideoReader(fileNames{file_i});
    %% Creates datetime timeseries for output file
    % Get end time of avi file
    video_dir = dir(fileNames{file_i});
    video_end_time = datetime([video_dir.date]);
    % Get the duration of the video in seconds
    videoLengthSeconds = vid.Duration;
    % Specify the frequency (1/15 of a second)
    frequency = seconds(1/15);
    % Calculate the number of time steps
    numTimeSteps = floor(videoLengthSeconds * 15); % 1/15 of a second
    % Create a time series of date times starting from (video end time - video duration) and ending at the video end time
    timeSeries = datetime(video_end_time - (0:numTimeSteps) * frequency, 'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
    timeSeries = flip(timeSeries);  % Flips time series so it ends at the end time

    nof = vid.NumberOfFrames;
    data = cell(nof, 1);
    coeff=1e-07; % conversion from pixel^2 to m^2
    %% Extracting time, mass, diameter, Area, density, evaporation time for indivisual snow/rain droplet
    parfor i = 1:nof
        frames = read(vid,i);
        img = rgb2gray(frames);
        img1=imcrop(img,[1 1 384 288]);% included kapton tape. Tweak parameters but same for ALTA exp.
         ref_temp1(i)=max(max(double(img1)));% plate temperature if hydrometeros fall on kapton tape
        img=imcrop(img,[1 27 384-1 288-27]); %ROI. Tweak parameter
        imgbw = img > 70;% removed below 70- thresholding -binary
        imgbw = imfill(imgbw, 'Holes');
         props = regionprops(imgbw, 'Centroid', 'Area','BoundingBox');
       if (isempty(props))
            continue;
        end
        bbobj = cat(1, props.BoundingBox);
        cobj = round(cat(1, props.Centroid)); % Centroids of snow
        aobj = cat(1, props.Area); % Areas of snow
        aobj=aobj* coeff;% conversion in m^2
        hei = bbobj(:, 3)*sqrt(coeff);% major axis
        wid = bbobj(:, 4)*sqrt(coeff);% minor axis
        ellip=hei.*wid; % ellipse area - Circumscribed area
        idx = sub2ind(size(img), cobj(:, 2), cobj(:, 1));
        iobj = double(img(idx));% Intensities of Centroid pixels of snow
       iobj=145-(iobj.*(145./ref_temp1(i)));%% 145 is max. temp range in screen setting
         prod=aobj.*iobj;% product of individual area and temp. diff.
        sa(i)=sum(aobj);% sum of area of hydrometeors in each frame
        s(i)=sum(prod);% sum of product of individual area and temp. diff in each frame
        cai = cat(2, cobj, aobj, iobj, ellip, hei,wid);
        data{i} = cai; % xc, yc, area, intensity
        % figure(2);
        % imshow(imgbw);
        % pause(0.0005);
        
    end
    %% SWE rate calculation
    A_hot=0.0101;
    fps=15;
    L=2.594*10^6;
    del=1/fps; 
    m1=0.0029*(s);% mass evaporates in each frame
    m1=m1/fps;
    swe=m1/A_hot;
    E1=m1*L./(sa);
    %% Sorting one frame to others frame
    data_sorted = cell(size(data));
    data_sorted{1} = data{1};
    
    for i = 2:nof
        data_sorted{i} = data{i};
        data_sorted{i} = sortPositions_v2(data_sorted{i-1}, data_sorted{i}, 20);
        
    end
    
    nmax = max(cellfun(@(x) size(x, 1), data_sorted, 'UniformOutput', 1));
    data_sorted = cellfun(@(x) cat(1, x, zeros(nmax-size(x, 1), 7)),...
        data_sorted, 'UniformOutput', 0);
    %% Testing
    
    % 3rd column is area and 4th column is intensity of centroid
    temparea = cellfun(@(x) x(:, 3), data_sorted, 'UniformOutput', 0);
    tempint = cellfun(@(x) x(:, 4), data_sorted, 'UniformOutput', 0);
    temparea1 = cellfun(@(x) x(:, 5), data_sorted, 'UniformOutput', 0);
    temphei = cellfun(@(x) x(:, 6), data_sorted, 'UniformOutput', 0);
    tempwid = cellfun(@(x) x(:, 7), data_sorted, 'UniformOutput', 0);
    droparea = cat(2, temparea{:});
    dropint = cat(2, tempint{:});
    droparea1 = cat(2, temparea1{:});
    drophei = cat(2, temphei{:});
    dropwid = cat(2, tempwid{:});
    % plot all area of snow between inial to final
    % figure, plot(transpose(droparea),'.-');
    
    droplife = {}; % Drop life
    dropmeanint1 = {}; % Drop centroid's mean intensity
    dropmeanint = {}; % Drop centroid's mean intensity
    dropmaxarea = {}; % Drop maximum area
    dropmaxarea1 = {};% Drop circumscribed area
    dropmaxhei = {};% Major axis
    dropmaxwid = {};% minor axis
    dropmass = {}; % Drop mass
    dropvolume = {}; % Drop volume
    dropdensity = {}; % Drop density
    time1=[];
    % time2=[];
    for i = 1:nmax
        
        temp = diff(droparea(i, :) > 0);
        idxbegin = find(temp > 0);
        idxend = find(temp < 0);% time when snow falls on the plate
        time1=cat(2,time1,idxbegin);% time when snow arrives the plate
        if isempty(idxend) % losing data at end those are not evaporating completely
            continue; 
        else
    %          time2=cat(2,time2,idxend);% time when snow arrives the plate
    %     if isempty(idxend) % losing data at end those are not evaporating completely
    %         continue; 
    %     else
            bw1temp = droparea(i, idxbegin(1):idxend(end)+1);
            bw2temp = bw1temp > 0;
            bw3temp = dropint(i, idxbegin(1):idxend(end)+1);
            bw4temp = droparea1(i, idxbegin(1):idxend(end)+1);
            bw5temp = drophei(i, idxbegin(1):idxend(end)+1);
            bw6temp = dropwid(i, idxbegin(1):idxend(end)+1);
            % anydrop which lives for less than 'tminlife' is discarded
            tminlife = 0;% no of frame
            bw2temp = bwareaopen(bw2temp, tminlife);
            % integration
            propstemp = regionprops(bw2temp, 'PixelIdxList');
        
            for j = 1:numel(propstemp)
                
                dropmaxarea{end+1} = max(bw1temp(propstemp(j).PixelIdxList));
                dropmaxarea1{end+1} = max(bw4temp(propstemp(j).PixelIdxList));
                dropmaxhei{end+1} = max(bw5temp(propstemp(j).PixelIdxList));
                dropmaxwid{end+1} = max(bw6temp(propstemp(j).PixelIdxList));
                droplife{end+1} = numel(propstemp(j).PixelIdxList);
                dropmeanint{end+1} = max(bw3temp(propstemp(j).PixelIdxList));
                dropmeanint1{end+1} = mean(bw3temp(propstemp(j).PixelIdxList));
                axitemp = bw1temp(propstemp(j).PixelIdxList).*...
                    bw3temp(propstemp(j).PixelIdxList);
                
                % the following trapz algorithm assumes that dt = 1
                %constant k/dLv = 0.0031
                dropmass{end+1} = trapz(0.0029*axitemp);
                % fps is not consider here
            end
        end
        
        end
    
    
    %% OUTPUT SWE .. calculation and save data
    % Drops present at the begining (t = 0) and ending (t = end) are excluded
    
    mass=cell2mat(dropmass);
    fps=15;%change fps here
    m=mass/fps;
    A=cell2mat(dropmaxarea);% actual area
    A1=cell2mat(dropmaxarea1);% ellipse area
    h1=cell2mat(dropmaxhei);% height
    w1=cell2mat(dropmaxwid);% Width
    D=1.12*A.^0.5; % diameter of drop
    D1=(6.*m/3140).^0.33;% water equi. diameter
    T=cell2mat(dropmeanint);% temperature diffrence between plate and water droplet
    T1=cell2mat(dropmeanint1);% temperature diffrence between plate and water droplet
    %T1=mean(T);
    t_evp=cell2mat(droplife)*(1/fps);% evaporation time
    tt=sort(time1')/fps;% time in sec
    % tt1=sort(time2')/fps;
    % tt=tt(1:length(m));% some snowflakes in between of evaporation- not counted
    V=0.75*A.^1.5;% Spherical volume
    rho_s=m./V; % density calculation: spherical assumption
    L=2.26*10^6;
    E=m*L./(A.*t_evp);%  heat flux method: energy per unit area per time
    d_heat1= 6.4418e+04*m./(A.*t_evp.*T);
     d_heat2= 6.4418e+04*m./(A.*t_evp.*T1);
    % rho_h=26*E/(5*10^4); %Density- heat flux method
    % rho_h1=E/729; %Density- heat flux method2
     v1=m./d_heat2;% volume of each snowflakes using heat flux method density
     hei=cell2mat(dropmaxhei);
     width=cell2mat(dropmaxwid);
    %% data all output
    jj= time1(1:length(m));
    Data1=[tt(1:length(m)) m' D' A' A1' t_evp' rho_s' d_heat1' d_heat2' E' v1' D1' hei' width' jj' T' T1'];
    Data2=[timeSeries' swe' E1'];
    % %  Data3=ref_temp1;

    %% Saves processed output data for single .avi file
    % Replace the file extension with '.txt'
    outputFileName = strrep(fileNames{file_i}, '.avi', '');
    disp(['Saved Output for: ', outputFileName])
    % Saves a file for time series and particle series
    save([outputDir, outputFileName, '_TS.txt'], 'Data2', '-ASCII'); 
    save([outputDir, outputFileName, '_Particles.txt'], 'Data1', '-ASCII');


     TotalSWE = sum(swe);
     density1=mean(d_heat1, 'omitnan');
     density2=mean(d_heat2, 'omitnan');
     TotalSnow1=1000* (TotalSWE./density1);
     TotalSnow2=1000* (TotalSWE./density2);
     
      %% Save Caclulated SWE and Density Data data and push to Eric's webpage
    
        % Calculate and store metrics for table:
        
        % Obtain date for each table row:
    %     Data_dir = dir('E:\DEID_Processor\Video_Compilations\Dec08_2023.avi');
    %     Data_Date = [Data_dir.date];
    %     date_time = datetime(Data_Date);
    %     
    %     format bank
    % 
    %     % Create a table with the stored data: 
    %     table_to_write = table(date_time, TotalSWE, density1, density2, TotalSnow1, TotalSnow2);
    % 
    %     % Rename the headers of the table:
    %     table_to_write.Properties.VariableNames = ["Date Time", "SWE(mm)", "Density1 (kg/m^3)", "Density2 (kg/m^3)", "SnowAccumulation1 (mm)", "SnowAccumulation2 (mm)"];
    % 
    %     % Write data to text file:
    %     writetable(table_to_write,'Dhiraj_DEID_Data.txt', 'Delimiter', '\t',  'WriteMode', 'Append');
     
      
    % end
    %    save('Feb_05_1_ref_temp.txt', 'Data3', '-ASCII')
    % % tt is time in second, m is mass in kg, D is eqv. circular diameter in m, 
    % A is actual area in m^2 , area1 is ellipse raea m^2, t_evp is evaporation
end