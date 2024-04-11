function sortPos = sortPositions_v2(Prev_Frame_Hydometer_Data, Cur_Frame_Hydometer_Data, dia)
% Prev_Frame_Hydometer_Data = The hydrometer data from the previous frame
% in the form Hydrometer_Centroid, Hydrometero_Area,
% [Plate_Hydrometer_DeltaT, Hydrometer_ellipse_area,Hydrometer_Major_Axis,
%Hydrometer_Minor_Axis]
% Cur_Frame_Hydometer_Data = The hydrometer data from the current frame
% in the form Hydrometer_Centroid, Hydrometero_Area,
% [Plate_Hydrometer_DeltaT, Hydrometer_ellipse_area,Hydrometer_Major_Axis,
%Hydrometer_Minor_Axis]
% dia = Some threshold for the change in the RMS between the centroid and
% area between succsive images
% DESSCRIPTION: This function... I'm not really sure yet...

temp1 = zeros(size(Cur_Frame_Hydometer_Data)); %temporary 1
temp2 = zeros(size(Cur_Frame_Hydometer_Data)); %temporary 2
non_zero_Previous_Frame = find(sum(Prev_Frame_Hydometer_Data, 2)); %Here, find just finds non-zero indices... Sum along all dim two? What the hell? then find?
non_zero_Current_Frame = find(sum(Cur_Frame_Hydometer_Data, 2));
N_Hydrometers_Prev_Frame = length(non_zero_Previous_Frame); % Number hydrometeros in previous non-zero frames 
N_Hydrometers_Cur_Frame = length(non_zero_Current_Frame); % Number hydrometeros in previous non-zero frames 

flag = false; %boolean flag..
kk = 1;

for ii = 1:N_Hydrometers_Cur_Frame % For each hydrometer in current frame.
    found = false;
    for jj = 1:N_Hydrometers_Prev_Frame % For each hydrometer in previous frame.
        %Here Dhiraj is looking for the rms for the change in the
        %[centroid, area] between succesive frames
        rms_hydrometer_delta_centroid_area = norm(Cur_Frame_Hydometer_Data(non_zero_Current_Frame(ii), 1:2)-...
            Prev_Frame_Hydometer_Data(non_zero_Previous_Frame(jj), 1:2)); 
        %if it changes less than threshold, the hydrometer is found?
        if (rms_hydrometer_delta_centroid_area <= dia)
            found = true; %set found to true
            temp1(non_zero_Previous_Frame(jj), :) = Cur_Frame_Hydometer_Data(non_zero_Current_Frame(ii), :);
            break %leave previous frame loop
        end
    end
    %no hydrometer found between frames
    if found == false
        flag = true;
        temp2(kk, :) = Cur_Frame_Hydometer_Data(non_zero_Current_Frame(ii), :);
        kk = kk+1;
    end
    
end

Cur_Frame_Hydometer_Data = temp1;

if flag == true
    idx = find(sum(Cur_Frame_Hydometer_Data, 2) == 0, kk-1);
    Cur_Frame_Hydometer_Data(idx, :) = temp2(1:kk-1, :);
end

sortPos = Cur_Frame_Hydometer_Data;