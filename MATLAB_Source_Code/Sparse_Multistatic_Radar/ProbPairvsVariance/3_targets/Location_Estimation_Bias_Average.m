
function [Estimated_Location_Average, Estimated_Location_Bias,True_Index_Record] = Location_Estimation_Bias_Average(Num_Targets, est_signal_1,est_signal_2, Target_Positions_pol, Target_Positions_car,pol_grid, x_vec, y_vec,I,J)

True_Index_Record = zeros(Num_Targets,2);
for index = 1:Num_Targets
    Tmp = sqrt((Target_Positions_pol(index,1) - pol_grid(:,1)).^2 + (Target_Positions_pol(index,2) - pol_grid(:,2)).^2);
    [~, Tmp_index_1] = min(Tmp);
    %Position_Record(index,:) = Target_Grids(Tmp_index,:);
    [True_Index_Record(index,2),True_Index_Record(index,1)] = ind2sub([I,J],Tmp_index_1);
end
% Bound is the search bound around each source in one side

N_Point = 5;


Estimated_Index_1  = zeros(Num_Targets,2);
Estimated_Index_2  = zeros(Num_Targets,2);

Estimated_Location_1  = zeros(Num_Targets,2);
Estimated_Location_2  = zeros(Num_Targets,2);
Estimated_Location_Average  = zeros(Num_Targets,2);
Estimated_Location_Bias = zeros(Num_Targets,1);
 %[Target_Positions_x,Target_Positions_y]=pol2car(Target_Positions_pol(1),Target_Positions_pol(2));
 for index = 1:Num_Targets
     %Use these when the search is focused on the estimated target
     %Index_Search_Range_r = (-N_Point:1:N_Point) +True_Index_Record(index,2);
     %Index_Search_Range_a = (-N_Point:1:N_Point) +True_Index_Record(index,2);
     Index_Search_Range_x = 1:1:length(x_vec);
     Index_Search_Range_y = (-N_Point:1:N_Point) +True_Index_Record(index,2);

     Location_Results_Search_1 = est_signal_1(Index_Search_Range_y,Index_Search_Range_x);
     Location_Results_Search_2 = est_signal_2(Index_Search_Range_y,Index_Search_Range_x);
     [~, Tmp_index_1] = max(Location_Results_Search_1(:));
     [~, Tmp_index_2] = max(Location_Results_Search_2(:));
     [Estimated_Index_1(index,2),Estimated_Index_1(index,1)] = ind2sub(size(Location_Results_Search_1),Tmp_index_1); %[range,doa]
     [Estimated_Index_2(index,2),Estimated_Index_2(index,1)] = ind2sub(size(Location_Results_Search_2),Tmp_index_2); %[range,doa]

     Estimated_Location_1(index,:) = [x_vec(Index_Search_Range_x(Estimated_Index_1(index,1))), y_vec(Index_Search_Range_y(Estimated_Index_1(index,2)))];
     Estimated_Location_2(index,:) = [x_vec(Index_Search_Range_x(Estimated_Index_2(index,1))), y_vec(Index_Search_Range_y(Estimated_Index_2(index,2)))];
     Estimated_Location_Average(index,:) = (Estimated_Location_1(index,:)+Estimated_Location_2(index,:))/2;
     
     Estimated_Location_Bias(index,:) = (Estimated_Location_Average(index,1)-Target_Positions_car(index,1)).^2+(Estimated_Location_Average(index,2)-Target_Positions_car(index,2)).^2;%sum of x and y bias
 end

end

