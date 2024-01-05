function [R,theta] = car2pol(x,y)
% This function converts Cartesian to polar coordinates using equation (24)
% Author: Ali Moussa - Last modified: 05/01/2024
R=sqrt(x.^2+y.^2);
theta=asind(x./R);  
end