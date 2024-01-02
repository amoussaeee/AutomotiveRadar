function [Estimated_Doppler,Doppler_Bias,True_Index_Record] = Doppler_Estimation_Bias(Num_Targets, est_signal, Target_Doppler,Doppler_vec)

True_Index_Record = zeros(Num_Targets,1);
for index = 1:Num_Targets
    Tmp = sqrt((Target_Doppler(index) - Doppler_vec).^2);
    [~, Tmp_index] = min(Tmp);
    True_Index_Record(index)=Tmp_index;
end

N_Point = 5;
Estimated_Index  = zeros(Num_Targets,1);
Estimated_Doppler  = zeros(Num_Targets,1);
%Actual_velocity=zeros(Num_Targets,1);
Doppler_Bias = zeros(Num_Targets,1);
for index = 1:Num_Targets
     %Use these when the search is focused on the estimated target
    Index_Search_Range= (-N_Point:1:N_Point) +True_Index_Record(index);
     Doppler_Results_Search = est_signal(Index_Search_Range);
     [~, Tmp_index] = max(Doppler_Results_Search(:));
     Estimated_Index(index)=Tmp_index;
     Estimated_Doppler(index) = Doppler_vec(Index_Search_Range(Estimated_Index(index)));
     %Actual_velocity(index)=Doppler2velocity(est_range_DOA(1),est_range_DOA(2),R_h,theta_h,Estimated_Doppler(index),v_s);
     Doppler_Bias(index) = (Estimated_Doppler(index)-Target_Doppler(index)).^2;
 end

end

