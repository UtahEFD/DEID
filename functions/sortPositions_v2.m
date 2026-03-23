function sortPos = sortPositions_v2(Prev_Frame_Hydometer_Data, Cur_Frame_Hydometer_Data, dia)
% **inputs**
% Prev_Frame_Hydometer_Data = hydrometer data from the previous frame
% Cur_Frame_Hydometer_Data = hydrometer data from the current frame
% hydrometer_data = (h_centroid, h_area, plate_h_dtemp, h_elipse_area, h_major_axis, h_minor_axis)
% dia = Some threshold for the change in the RMS between the centroid and
% area between succsive images
% **description**
% match hydrometeors between two consecutive frames of an image sequence, based on their centroid positions.
% if a hydrometeor in the current frame can be matched to one in the previous frame, it is kept in the same row position.
% if a hydrometeor cannot be matched (unique centroid), it gets stored separately and re-inserted into current frame array at the first empty rows.
% i.e the rows of Cur_Frame_Hydometer_Data correspond as much as possible to the same physical hydrometeor across frames.

temp1 = zeros(size(Cur_Frame_Hydometer_Data)); %temporary 1
temp2 = zeros(size(Cur_Frame_Hydometer_Data)); %temporary 2
non_zero_Previous_Frame = find(sum(Prev_Frame_Hydometer_Data, 2)); %Here, find just finds non-zero indices... Sum along all dim two? What the hell? then find?
non_zero_Current_Frame = find(sum(Cur_Frame_Hydometer_Data, 2));
N_Hydrometers_Prev_Frame = length(non_zero_Previous_Frame); % Number hydrometeros in previous non-zero frames 
N_Hydrometers_Cur_Frame = length(non_zero_Current_Frame); % Number hydrometeros in previous non-zero frames 

flag = false; %boolean flag..
kk = 1;

for ii = 1:N_Hydrometers_Cur_Frame % for each hydrometer in current frame.
    found = false;
    for jj = 1:N_Hydrometers_Prev_Frame % for each hydrometer in previous frame.
        % looks for rms of change in the centroid between succesive frames
        rms_hydrometer_delta_centroid_area = norm(Cur_Frame_Hydometer_Data(non_zero_Current_Frame(ii), 1:2)-...
            Prev_Frame_Hydometer_Data(non_zero_Previous_Frame(jj), 1:2)); 
        % if centroid difference is below threshold (dia), consider it the
        % same hydrometor and store in temp1; break loop 
        if (rms_hydrometer_delta_centroid_area <= dia)
            found = true;
            temp1(non_zero_Previous_Frame(jj), :) = Cur_Frame_Hydometer_Data(non_zero_Current_Frame(ii), :);
            break 
        end
    end
    % if centroid difference is above threshold, consider it new
    % hydrometeor and store in temp2
    if found == false
        flag = true;
        temp2(kk, :) = Cur_Frame_Hydometer_Data(non_zero_Current_Frame(ii), :);
        kk = kk+1;
    end
    
end

Cur_Frame_Hydometer_Data = temp1; % starting with only matched hydrometeors 
% if there were unmatched (new) hydrometeors, flag=true, inserts into empty
% rows of cur_frame_hydrometeor_data

if flag == true
    idx = find(sum(Cur_Frame_Hydometer_Data, 2) == 0, kk-1);
    Cur_Frame_Hydometer_Data(idx, :) = temp2(1:kk-1, :);
end

sortPos = Cur_Frame_Hydometer_Data;