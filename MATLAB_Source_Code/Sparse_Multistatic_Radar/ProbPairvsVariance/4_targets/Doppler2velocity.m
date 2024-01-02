function [V_k] = Doppler2velocity(R_k,theta_k,R_s,theta_s,V_k,v_s)
%This function converts the Doppler of the vehicle into the velocity term
%in equation 5 for v_k.
gamma_k=abs(theta_s-theta_k);
hat_R_k=sqrt(R_k.^2+R_s^2-2.*R_k*R_s.*cosd(abs(theta_s-theta_k)))+R_k;
alpha_k=asind(R_s.*sind(gamma_k)./(hat_R_k-R_k))+theta_k;
V_k=(V_k+v_s.*cosd(theta_k))./(cosd(alpha_k)+cosd(theta_k));

end