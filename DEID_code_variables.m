clearvars;
close all
%% Extracting mass, diameter, Area, density, evaporation time for indivisual snow/rain droplet 
nof = 3300; % Number of frames
data = cell(nof, 1);
% Subtract background to remove noise in the hotplate
% background image (any image may ne first) without hydrometeors
% imgbg = imread(['~/Downloads/40ul/img_', num2str(0, '%04d'), '.tif']);

for i = 2:nof

    imgname = ['/Users/apple/Documents/TARP_DEID/Jan_11_good_snow_/fps_12_03_0001_',...
                        num2str(i-1, '%04d'), '.jpg'];
    % Subtracting the background image is not necessary
%      img = imread(imgname)-imgbg;
    img = rgb2gray(imread(imgname));
    img1=imcrop(img,[2 2 463 357]);% included kapton tape
%     ref_temp(i)=mean(double(img1));
    ref_temp1(i)=max(max(double(img1)));% if hydrometeros fall on kapton tape
    % and Emissivity role snow or water > than kapton tape 
    img=imcrop(img,[1 47 463-1 357-47]); %ROI
    imgbw = img > 70;% removed below 70
    imgbw = imfill(imgbw, 'Holes');
    imgbw = bwareaopen(imgbw, 10);
    imgbw = imclearborder(imgbw);
    props = regionprops(imgbw, 'Centroid', 'Area');
    if (isempty(props))
        continue;
    end
    cobj = round(cat(1, props.Centroid)); % Centroids of objects
    aobj = cat(1, props.Area); % Areas of objects
    aobj=aobj* 1e-07;
    idx = sub2ind(size(img), cobj(:, 2), cobj(:, 1));
    iobj = double(img(idx));% Intensities of Centroid pixels of objects
    iobj=(ref_temp1(i)-iobj)*(105/255);% diffrence between plate and water top
    cai = cat(2, cobj, aobj, iobj);
    data{i} = cai; % xc, yc, area, intensity
    % figure(2); 
    % imshow(imgbw);
    % pause(0.5);
    
end

%% Sorting one frmae to others frame
data_sorted = cell(size(data));
data_sorted{1} = data{1};

for i = 2:nof
    
    data_sorted{i} = data{i};
    data_sorted{i} = sortPositions_v2(data_sorted{i-1}, data_sorted{i}, 20);
    
end

nmax = max(cellfun(@(x) size(x, 1), data_sorted, 'UniformOutput', 1));
data_sorted = cellfun(@(x) cat(1, x, zeros(nmax-size(x, 1), 4)),...
                        data_sorted, 'UniformOutput', 0);
%% Testing

% 3rd column is area and 4th column is intensity of centroid
temparea = cellfun(@(x) x(:, 3), data_sorted, 'UniformOutput', 0);
tempint = cellfun(@(x) x(:, 4), data_sorted, 'UniformOutput', 0);
droparea = cat(2, temparea{:});
dropint = cat(2, tempint{:});
% plot all area of snow between inial to final
plot(transpose(droparea),'.-');

droplife = {}; % Drop life
dropmeanint = {}; % Drop centroid's mean intensity
dropmaxarea = {}; % Drop maximum area
dropmass = {}; % Drop mass
dropvolume = {}; % Drop volume
dropdensity = {}; % Drop density


for i = 1:nmax
    
    temp = diff(droparea(i, :) > 0);
    idxbegin = find(temp > 0);
    idxend = find(temp < 0); 
    
    bw1temp = droparea(i, idxbegin(1):idxend(end)+1);
    bw2temp = bw1temp > 0;
    bw3temp = dropint(i, idxbegin(1):idxend(end)+1);
    
    % anydrop which lives for less than 'tminlife' is discarded
    tminlife = 3;
    bw2temp = bwareaopen(bw2temp, tminlife);
    % integration
    propstemp = regionprops(bw2temp, 'PixelIdxList');
        
        for j = 1:numel(propstemp)
            
            dropmaxarea{end+1} = max(bw1temp(propstemp(j).PixelIdxList));
            droplife{end+1} = numel(propstemp(j).PixelIdxList);
            dropmeanint{end+1} = mean(bw3temp(propstemp(j).PixelIdxList));
            
            axitemp = bw1temp(propstemp(j).PixelIdxList).*...
                        bw3temp(propstemp(j).PixelIdxList);
                    
            % the following trapz algorithm assumes that dt = 1
            %constant k/dLv = 
            dropmass{end+1} = trapz(0.0031*axitemp);
            % fps is not consider here
        end
    
end

dropvolume = cellfun(@(x,y,z) x*y*z, droplife, dropmaxarea, dropmeanint, ...
                        'UniformOutput', 0);
dropdensity = cellfun(@(x,y) x/y, dropmass, dropvolume, ...
                        'UniformOutput', 0);
%% OUTPUT SWE .. calculation and save data
Drops present at the begining (t = 0) and ending (t = end) are excluded
%% output Run section after first one
mass=cell2mat(dropmass);
% for 20 fps time in each frame = 1/20 sec.
fps=20;%change fps here
m=mass/fps;
%area
A=cell2mat(dropmaxarea);
%D_eff
D=1.12*A.^0.5;
% D_eff for rain droplets
D1=(6.8.*m/3140).^0.33;
% temperature diffrence between plate and water droplet
T=cell2mat(dropmeanint);
T1=mean(T);
% evaporation time
t=cell2mat(droplife)*(1/fps);
% volume
V=0.75*A.^1.5;
% density calculation: spherical assumption
rho_s=m./V;
L=2.26*10^6;
% heat flux method: energy per unit area,time
E=m*L./(A.*t);
%Linear relation between heat flux and density
% need more vaification: density heat flux method
 rho_h=26*E/(5*10^4);
 %% SWE calculation
 % area of hot plate in m^2: calibration parameter
A_hot=size(img,1)*size(img,2)*1e-07;
% choose the duration of SWE
nof=3300;
time=nof./fps;
time=time/3600;
m1= sum(m);
% SWE rate in mm/hr
SWE=m1./(time.*A_hot);
% snow rate in MM/hr
rate_snow=SWE*1000/mean(rho_s);
% water accmmulation
Accum_SWE=SWE*time;
% snow accumulation
Accum_snow=rate_snow*time;
% saving data: all input
Data1=[m' D' A' t' rho_s' rho_h' E'];
Data2=[SWE rate_snow Accum_SWE Accum_snow T1];
figure, plot(m,'o')
figure,plot(rho_h,'o')
figure,plot(rho_s,'o')

save('mass_Dia_Area_evptime_sphdensity_heatdensity_DelT_rain_001.txt', 'Data1','-ASCII')

  save('SWE_snowrate_SWEaccu_snowaccu_00510.txt', 'Data2','-ASCII')

% write all information that need to re-plot
% from the analysis because of the incomplete information or process.
%fps 12, T_plate = 85C, k/(dLv) = 4.1/2.670, pix2m=1.99e-4, crop = [34 62 327-34 268-62]
