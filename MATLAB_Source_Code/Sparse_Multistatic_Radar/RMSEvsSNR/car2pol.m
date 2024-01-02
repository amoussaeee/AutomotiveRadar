function [R,theta] = car2pol(x,y)
R=sqrt(x.^2+y.^2);
theta=asind(x./R);  
end