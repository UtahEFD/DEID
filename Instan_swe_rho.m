%% This code for instantaneous SWE and average density calculation
clear all
clearvars;
%% Extracting positions 
nof = 3300; % Number of frames
data = cell(nof, 1);
% background image (first image)
% imgbg = imread(['~/Downloads/40ul/img_', num2str(0, '%04d'), '.tif']);

for i = 2:nof

    imgname = ['/Users/apple/Documents/TARP_DEID/New folder (4)/AAAAA/fps_12_04_1630_0070_',...
                        num2str(i-1, '%04d'), '.jpg'];
    % Subtracting the background image is not necessary: only when static
    % mark on the plate
    % img = imread(imgname)-imgbg;
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
    prod=aobj.*iobj;% product of individual area and temp. diff.
    sa(i)=sum(aobj);% sum of area of hydrometeors in each frame
    s(i)=sum(prod);% sum of product of individual area and temp. diff in each frame
    cai = cat(2, cobj, aobj, iobj);
    data{i} = cai; % xc, yc, area, intensity
    % figure(2); 
    % imshow(imgbw);
    % pause(0.5);
    
end
%% Iinstantaneous SWE calculation
L=2.26*10^6;% latent heat of vaparization
fps=12;
del=1/fps;
t=0:del:del*3299;% time in sec
m=0.0031*s;% mass evaporates in each frame
E=m*L./(sa);% Heat flux in each frame
rho=45*E./(5*10^4);
A_hot=size(img,1)*size(img,2)*1e-07;% plate area in m^2
swe=3600*m./A_hot; %mm/hr instantaneous swe rate 
t=t/3600;% time in hr
total_acc = trapz(t, swe);% water accumulation in t (275 sec)time.
%% swe avarging over 10 sec period

 rho_10sec = zeros(27,1);
 for i=1:27
     rho_10sec(i,:)=mean(rho(120*(i-1)+1:120*(i-1)+120));
 end
     figure,  plot(rho_10sec,'o-')
%     figure,plot(t,E,'.-')

%% data saving
data_swe=[swe E];
save('mass_Dia_Area_evptime_sphdensity_heatdensity_DelT_rain_001.txt', 'data_swe','-ASCII')


