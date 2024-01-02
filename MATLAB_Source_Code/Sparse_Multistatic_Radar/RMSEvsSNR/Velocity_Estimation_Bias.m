function [Estimated_Velocity,Velocity_Bias,True_Index_Record] = Velocity_Estimation_Bias(Num_Targets, est_signal, Target_Velocity,velocity_vec)

True_Index_Record = zeros(Num_Targets,1);
for index = 1:Num_Targets
    Tmp = sqrt((Target_Velocity(index) - velocity_vec).^2);
    [~, Tmp_index] = min(Tmp);
    True_Index_Record(index)=Tmp_index;
end

N_Point = 5;
Estimated_Index  = zeros(Num_Targets,1);
Estimated_Velocity  = zeros(Num_Targets,1);
%Actual_velocity=zeros(Num_Targets,1);
Velocity_Bias = zeros(Num_Targets,1);
for index = 1:Num_Targets
     %Use these when the search is focused on the estimated target
    Index_Search_Range= [(-N_Point:1:N_Point) +True_Index_Record(index);
     Velocity_Results_Search = est_signal(Index_Search_Range);
     [~, Tmp_index] = max(Velocity_Results_Search(:));
     Estimated_Index(index)=Tmp_index;
     Estimated_Velocity(index) = velocity_vec(Index_Search_Range(Estimated_Index(index)));
     %Actual_velocity(index)=Doppler2velocity(est_range_DOA(1),est_range_DOA(2),R_h,theta_h,Estimated_Doppler(index),v_s);
     Velocity_Bias(index) = (Estimated_Velocity(index)-Target_Velocity(index)).^2;
 end

end

