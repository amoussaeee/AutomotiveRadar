function [Actual_Doppler, Doppler_Bias,Index_Record] = doppler_estimation_errors(Num_Targets, est_signal, Target_Doppler,est_range_DOA, Doppler_vec,G_d,R_s,theta_s,v_s)

Index_Record = zeros(Num_Targets,1);
for index = 1:Num_Targets
    Tmp = sqrt(Target_Doppler(index) - Doppler_vec).^2 ;
    [~, Tmp_index] = min(Tmp);
    Index_Record(index)=Tmp_index;
end


Actual_Index  = zeros(Num_Targets,1);
Actual_Doppler  = zeros(Num_Targets,1);
Actual_velocity=zeros(Num_Targets,1);
Doppler_Bias = zeros(Num_Targets,1);
for index = 1:Num_Targets
     Index_Search_Range = 1:1:length(Doppler_vec);
     Doppler_Results_Search = est_signal(Index_Search_Range);
     [~, Tmp_index] = max(Doppler_Results_Search(:));
     Actual_Index(index)=Tmp_index;
     Actual_Doppler(index) = Doppler_vec(Index_Search_Range(Actual_Index(index)));
     Actual_velocity(index)=doppler2velocity(est_range_DOA(1),est_range_DOA(2),R_s,theta_s,Actual_Doppler(index),v_s);
     Doppler_Bias(index) = (Actual_Doppler(index)-Target_Doppler(index)).^2;
 end

end

