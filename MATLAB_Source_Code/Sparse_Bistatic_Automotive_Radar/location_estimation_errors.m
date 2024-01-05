function [Actual_Location, Location_Bias,Index_Record] = location_estimation_errors(Num_Targets, est_location, Target_Positions_pol ,pol_grid, x_vec, y_vec,I,J)

Index_Record = zeros(Num_Targets,2);
for index = 1:Num_Targets
    Tmp = sqrt((Target_Positions_pol(index,1) - pol_grid(:,1)).^2 + (Target_Positions_pol(index,2) - pol_grid(:,2)).^2);
    [~, Tmp_index] = min(Tmp);
    [Index_Record(index,2),Index_Record(index,1)] = ind2sub([I,J],Tmp_index);
end

Actual_Index  = zeros(Num_Targets,2);
Actual_Location  = zeros(Num_Targets,2);
Location_Bias = zeros(Num_Targets,1);
[Target_Positions_x,Target_Positions_y]=pol2car(Target_Positions_pol(1),Target_Positions_pol(2));
 for index = 1:Num_Targets
     Index_Search_Range_x = 1:1:length(x_vec);
     Index_Search_Range_y = 1:1:length(y_vec);
     Location_Results_Search = est_location(Index_Search_Range_y,Index_Search_Range_x);
     [~, Tmp_index] = max(Location_Results_Search(:));
     [Actual_Index(index,2),Actual_Index(index,1)] = ind2sub(size(Location_Results_Search),Tmp_index); %[range,doa]
     Actual_Location(index,:) = [x_vec(Index_Search_Range_x(Actual_Index(index,1))), y_vec(Index_Search_Range_y(Actual_Index(index,2)))];
     Location_Bias(index,:) = (Actual_Location(index,1)-Target_Positions_x(index)).^2+(Actual_Location(index,2)-Target_Positions_y(index)).^2;%sum of x and y bias
 end

end

