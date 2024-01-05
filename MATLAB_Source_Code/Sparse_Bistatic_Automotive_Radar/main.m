%% Introduction
% This program is to simulate an automotive bistatic scenario where the
% 2D location is estimated followed by Doppler estimation using l21 
% minimisation under the a group sparsity framework. A more detailed
% description of the method used here is published in "A Two-Stage 
% Sparsity-Based Method for Location and Doppler Estimation in Bistatic 
% Automotive Radar".

% Author: Ali Moussa - Last modified: 05/01/2024
%% Initialisation
clear all;
close all
format long;
format compact;
tic
fprintf('Program executed at:')
fprintf(string(datetime))
%rng(20,'philox')

%%% Define Chirp Parameters %%%
fc  = 77e09;          %centre frequency
B   = 150e6;          %chirp bandwidth 
Tc  = 30e-6;          %modulation period 
c   = 3e08;           %speed of light
mu  = B/Tc;           %modulation rate
K   = 1;              %number of targets
Q   = [1 4 8];        %number of processed columns in group sparse

%%% Define  Clock Settings %%%
fs  = 5e6;            %sampling frequency (~bandwidth in complex sampling)
Ts  = 1/fs;           %sampling time
N   = floor(Tc/Ts);   %number of samples in each pulse

%%% Maximum Unambiguous Range %%%
Rm  = c*fs/mu;        %max detectable range (unused)

%%% Define ULA Parameters %%%
lambda  = c/fc;       %wavelength
d       = lambda/2;   %antenna spacing
L       = 9;          %number of antennas

%%% Define Slow Time Parameters %%%
M   = 128;         %number of chirps   
T   = 35e-6;       %pulse repetition period
v_min   = 20; %m/s
v_max   = 35; %m/s

%%% Searching Vehicle and Roadside Sensor Init %%%
R_s     = 30;       %metres
theta_s = -7.662;   %degrees
v_s     = 25;       %m/s


%%% 2D Rectangular Search Grid %%%
I       = 30;
J       = 7;
x_vec   = linspace(0,6,J);
y_vec   = linspace(74,45,I);
step_x  = abs(x_vec(2)-x_vec(1));
step_y  = abs(y_vec(2)-y_vec(1));

rect_grid_x = zeros(I,J);
rect_grid_y = zeros(I,J);

for j=1:J
    rect_grid_x(:,j)=x_vec(j);
end
for i=1:I
    rect_grid_y(i,:)=y_vec(i);
end

%%% Convert 2D rectangular grid to 1D polar grid of length G_p=IxJ %%%
G_p             = I*J;
pol_grid_R      = zeros(1,G_p);
pol_grid_theta  = zeros(1,G_p);
buffer_x        = rect_grid_x(:);
buffer_y        = rect_grid_y(:);

for g=1:G_p
    pol_grid_R(g)=sqrt(buffer_x(g)^2+buffer_y(g)^2);
    pol_grid_theta(g)=asind(buffer_x(g)/sqrt(buffer_x(g)^2+buffer_y(g)^2));%Automatically 0 deg when x=0 m.
end
pol_grid=[pol_grid_R',pol_grid_theta'];

%%% Define received signal power and noise power %%%
EIRP    = 1.995262315;  %Watts (eq to 33dBm)
Gr      = 10^(16/10);   %linear
RCS     = 1;            %linear (eq to 0dBsm)

%% Signal Model
q_k = @(u,a) sqrt((EIRP*Gr*(lambda^2)*RCS)/((4*pi)^3*(u^2)*(sqrt(u^2+R_s^2-2*u*R_s*cosd(abs(theta_s-a)))^2)))...
      *exp(-1i*2*pi*fc*(sqrt(u^2+R_s^2-2*u*R_s*cosd(abs(theta_s-a)))+u)/c);
p_k = @(u,a) kron(exp((-1i*2*pi*mu*(sqrt(u^2+R_s^2-2*u*R_s*cosd(abs(theta_s-a)))+u)*(0:N-1))/(c*fs)).',exp(-1i*2*pi*fc*d*sind(a)*(0:L-1)/c).'); 
v_k = @(v) exp(-1i*2*pi*fc*v*(0:M-1)*T/c);

%%% Constructing a LNxG_p location steering matrix %%%
P_g         = zeros(L*N,G_p);
V_min_var   = [];
V_max_var   = [];
for g=1:G_p
  P_g(:,g)  = p_k(pol_grid_R(g),pol_grid_theta(g));
  V_min_var = [V_min_var;velocity2doppler(pol_grid_R(g),pol_grid_theta(g),R_s,theta_s,v_min,v_s)];
  V_max_var = [V_max_var;velocity2doppler(pol_grid_R(g),pol_grid_theta(g),R_s,theta_s,v_max,v_s)];
end

%%% Constructing a MxG_d Doppler steering matrix %%%
G_d     = 64;
V_min   = min(V_min_var);
V_max   = max(V_max_var);
Doppler_vec  = linspace(V_min,V_max,G_d);
step_Doppler = abs(Doppler_vec(2)-Doppler_vec(1));
V_g=zeros(M,G_d);
for g=1:G_d
  V_g(:,g)=v_k(Doppler_vec(g));
end

%% SNR and Reconstruction Errors
%%% SNR %%%
SNR_vector  = 128:4:144; %dB
SNR_length=length(SNR_vector);

%%% Reconstruction Errors %%%
eb_vector_loc_1      = [0.0000211 0.0000131 0.0000084 0.0000053 0.0000034];
eb_vector_Doppler_1  = [0.0000060 0.0000045 0.0000035 0.0000035 0.0000030];

eb_vector_loc_4      = [0.0000408 0.0000253 0.0000159 0.0000100 0.0000064];
eb_vector_Doppler_4  = [0.0000120 0.0000065 0.0000070 0.0000055 0.0000045];

eb_vector_loc_8      = [0.0000564 0.0000356 0.0000225 0.0000141 0.0000089];
eb_vector_Doppler_8  = [0.0000170 0.0000090 0.0000090 0.0000070 0.0000060];


%% Loops
Loop_Num    = 1; %number of Monte Carlo Trials
Methods_num = 4;

%%% Buffer for mean squared error %%%
MSE_loc      = zeros(Methods_num,SNR_length);
MSE_Doppler  = zeros(Methods_num,SNR_length);

%%% Buffer for Sparse method (1 snapshot) %%%
loc_Bias_Store_CS_1         = zeros(K*Loop_Num,SNR_length);
Doppler_Bias_Store_CS_1     = zeros(K*Loop_Num,SNR_length);

%%% Buffer for Sparse method (4 snapshots) %%%
loc_Bias_Store_CS_4         = zeros(K*Loop_Num,SNR_length);
Doppler_Bias_Store_CS_4     = zeros(K*Loop_Num,SNR_length);

%%% Buffer for Sparse method (8 snapshots) %%%
loc_Bias_Store_CS_8         = zeros(K*Loop_Num,SNR_length);
Doppler_Bias_Store_CS_8     = zeros(K*Loop_Num,SNR_length);

%%% Buffer for MUSIC method (8 snapshots) %%%
loc_Bias_Store_SS_8         = zeros(K*Loop_Num,SNR_length);
Doppler_Bias_Store_SS_8     = zeros(K*Loop_Num,SNR_length);

%%% Buffer to investigate CVX Solver Record %%%
solver_record_loc       = [];
solver_record_Doppler   = [];



%% Monte Carlo Loop Start
for loop_index = 1:Loop_Num

fprintf('Current Monte Carlo Trial Number:')
fprintf(string(loop_index))

%% Target Parameters
%%% Draw target Cartesian coordinates from a uniform distribution %%%
x_k1 = unifrnd(x_vec(4)+step_x/4,x_vec(5)-step_x/4);        
y_k1 = unifrnd(y_vec(16)+step_y/4,y_vec(15)-step_y/4);      

%%% Convert target Cartesian coordinates to polar coordinates %%%
[R_k1,theta_k1] = car2pol(x_k1,y_k1);

%%% Draw target Doppler from a uniform distribution %%%
V_k1 = unifrnd(Doppler_vec(36)+step_Doppler/4,Doppler_vec(37)-step_Doppler/4);    

%%% Convert target Doppler to forward velocity %%%
v_k1 = doppler2velocity(R_k1,theta_k1,R_s,theta_s,V_k1,v_s);

Target_Position = [R_k1 theta_k1];
Target_Doppler  = V_k1;
Num_Targets     = K;

%% Generate a Noisy Raw Signal 
%%% Raw Signal %%%
Z_1= q_k(R_k1,theta_k1) * p_k(R_k1,theta_k1) * v_k(V_k1);

%%% Noise Init %%%
W_1 = randn([L*N M]);

%% SNR Loop
for SNR_index = 1:SNR_length
    SNR = SNR_vector(SNR_index);
    eb_loc_1      = eb_vector_loc_1(:,SNR_index);
    eb_Doppler_1  = eb_vector_Doppler_1(:,SNR_index);
 
    eb_loc_4      = eb_vector_loc_4(:,SNR_index);
    eb_Doppler_4  = eb_vector_Doppler_4(:,SNR_index);
 
    eb_loc_8      = eb_vector_loc_8(:,SNR_index);
    eb_Doppler_8  = eb_vector_Doppler_8(:,SNR_index);
    A_Noise = sqrt(EIRP^2/(10^(SNR/10))/2);

%%% Noise %%%
W_1_n = A_Noise*W_1;
W_1_t = W_1_n.'; %transpose of noise signal

%%% Noisy Signal %%%
Z_1_n = real(Z_1)+W_1_n;

%% Test values for adjusting the reconstruction error
SNR_test = snr(real(Z_1),real(W_1_n)); %unused

noise_test_loc_1     = norm(W_1_n(:,Q(1)),'fro');
noise_test_loc_4     = norm(W_1_n(:,Q(2)),'fro');
noise_test_loc_8     = norm(W_1_n(:,Q(3)),'fro');

noise_test_Doppler_1 = norm(W_1_t(:,Q(1)),'fro');
noise_test_Doppler_4 = norm(W_1_t(:,Q(2)),'fro');
noise_test_Doppler_8 = norm(W_1_t(:,Q(3)),'fro');

%% Q=1 %%
%%% Solve the location optimisation problem using l21 minimisation %%%
[est_location_signal,solver_status_loc]=l21_minimisation(Z_1_n,P_g,eb_loc_1,G_p,Q(1));
buffer_loc=reshape(est_location_signal,I,J);
if strcmp(solver_status_loc,'Solved')==0
   solver_record_loc=[solver_record_loc; [solver_status_loc string(SNR) string(eb_loc_1) string(noise_test_loc_1) string(loop_index) string(Q(1))]];
end
[est_loc_Value, loc_Bias,~] = location_estimation_errors(Num_Targets, buffer_loc, Target_Position, pol_grid, x_vec, y_vec,I,J); %Bias is error squared
loc_Bias_Store_CS_1(loop_index,SNR_index)=loc_Bias;
[est_R_k, est_theta_k]=car2pol(est_loc_Value(1),est_loc_Value(2));

%%% Plotting %%%
% Est_Surf_loc=buffer_loc;
% figure(1)
% imagesc(x_vec,y_vec,Est_Surf_loc)
% xlabel('x');
% ylabel('y');

%%% Solve the Doppler optimisation problem using l21 minimisation %%%
[est_Doppler_signal,solver_status_Doppler]=l21_minimisation(Z_1_n.',V_g,eb_Doppler_1,G_d,Q(1));
if strcmp(solver_status_Doppler,'Solved')==0
 solver_record_Doppler=[solver_record_Doppler; [solver_status_Doppler string(SNR) string(eb_Doppler_1) string(noise_test_Doppler_1) string(loop_index) string(Q(1))]];
end

[~, Doppler_Bias,~] = doppler_estimation_errors(Num_Targets, est_Doppler_signal, Target_Doppler,[est_R_k est_theta_k], Doppler_vec,G_d,R_s,theta_s,v_s); %Bias is error squared
Doppler_Bias_Store_CS_1(loop_index,SNR_index)=Doppler_Bias;

%%% Plotting %%%
% Est_Surf_Doppler=est_Doppler_signal;
% figure(2)
% plot(Doppler_vec,Est_Surf_Doppler)
% xlabel('Doppler');
% ylabel('Signal Level');


%% Q=4 %%
%%% Solve the location optimisation problem using l21 minimisation %%%
[est_location_signal,solver_status_loc]=l21_minimisation(Z_1_n,P_g,eb_loc_4,G_p,Q(2));
buffer_loc=reshape(est_location_signal,I,J);
if strcmp(solver_status_loc,'Solved')==0
 solver_record_loc=[solver_record_loc; [solver_status_loc string(SNR) string(eb_loc_4) string(noise_test_loc_4) string(loop_index) string(Q(2))]];
end
[est_loc_Value, loc_Bias,~] = location_estimation_errors(Num_Targets, buffer_loc, Target_Position, pol_grid, x_vec, y_vec,I,J); %Bias is error squared
loc_Bias_Store_CS_4(loop_index,SNR_index)=loc_Bias;
[est_R_k, est_theta_k]=car2pol(est_loc_Value(1),est_loc_Value(2));

%%% Plotting %%%
% Est_Surf_loc=buffer_loc;
% figure(3)
% imagesc(x_vec,y_vec,Est_Surf_loc)
% xlabel('x');
% ylabel('y');

%%% Solve the Doppler optimisation problem using l21 minimisation %%%
[est_Doppler_signal,solver_status_Doppler]=l21_minimisation(Z_1_n.',V_g,eb_Doppler_4,G_d,Q(2));
if strcmp(solver_status_Doppler,'Solved')==0
 solver_record_Doppler=[solver_record_Doppler; [solver_status_Doppler string(SNR) string(eb_Doppler_4) string(noise_test_Doppler_4) string(loop_index) string(Q(2))]]; %#ok<*AGROW> 
end
[~, Doppler_Bias,~] = doppler_estimation_errors(Num_Targets, est_Doppler_signal, Target_Doppler,[est_R_k est_theta_k], Doppler_vec,G_d,R_s,theta_s,v_s); %Bias is error squared
Doppler_Bias_Store_CS_4(loop_index,SNR_index)=Doppler_Bias;

%%% Plotting %%%
% Est_Surf_Doppler=est_Doppler_signal;
% figure(4)
% plot(Doppler_vec,Est_Surf_Doppler)
% xlabel('Doppler');
% ylabel('Signal Level');

%% Q=8 %%
%%% Solve the location optimisation problem using l21 minimisation %%%
[est_location_signal,solver_status_loc]=l21_minimisation(Z_1_n,P_g,eb_loc_8,G_p,Q(3));
buffer_loc=reshape(est_location_signal,I,J);
if strcmp(solver_status_loc,'Solved')==0
 solver_record_loc=[solver_record_loc; [solver_status_loc string(SNR) string(eb_loc_8) string(noise_test_loc_8) string(loop_index) string(Q(3))]];
end
[est_loc_Value, loc_Bias,~] = location_estimation_errors(Num_Targets, buffer_loc, Target_Position, pol_grid, x_vec, y_vec,I,J); %Bias is error squared
loc_Bias_Store_CS_8(loop_index,SNR_index)=loc_Bias;
[est_R_k, est_theta_k]=car2pol(est_loc_Value(1),est_loc_Value(2));

%%% Plotting %%%
% Est_Surf_loc=buffer_loc;
% figure(5)
% imagesc(x_vec,y_vec,Est_Surf_loc)
% xlabel('x');
% ylabel('y');

%%% Solve the Doppler optimisation problem using l21 minimisation %%%
[est_Doppler_signal,solver_status_Doppler]=l21_minimisation(Z_1_n.',V_g,eb_Doppler_8,G_d,Q(3));
if strcmp(solver_status_Doppler,'Solved')==0
 solver_record_Doppler=[solver_record_Doppler; [solver_status_Doppler string(SNR) string(eb_Doppler_8) string(noise_test_Doppler_8) string(loop_index) string(Q(3))]];
end
[~, Doppler_Bias,~] = doppler_estimation_errors(Num_Targets, est_Doppler_signal, Target_Doppler,[est_R_k est_theta_k], Doppler_vec,G_d,R_s,theta_s,v_s); %Bias is error squared
Doppler_Bias_Store_CS_8(loop_index,SNR_index)=Doppler_Bias;

%%% Plotting %%%
% Est_Surf_Doppler=est_Doppler_signal;
% figure(6)
% plot(Doppler_vec,Est_Surf_Doppler)
% xlabel('Doppler');
% ylabel('Signal Level');

%% MUSIC (Q=8) %%
est_location_signal=MUSIC(Z_1_n,real(P_g),Num_Targets,G_p,Q(3));
buffer_loc=reshape(est_location_signal,I,J);
[est_loc_Value, loc_Bias,~] = location_estimation_errors(Num_Targets, buffer_loc, Target_Position, pol_grid, x_vec, y_vec,I,J); %Bias is error squared
loc_Bias_Store_SS_8(loop_index,SNR_index)=loc_Bias;
[est_R_k, est_theta_k]=car2pol(est_loc_Value(1),est_loc_Value(2));

est_Doppler_signal=MUSIC(Z_1_n.',real(V_g),Num_Targets,G_d,Q(3));
[~, Doppler_Bias,~] = doppler_estimation_errors(Num_Targets, est_Doppler_signal, Target_Doppler,[est_R_k est_theta_k], Doppler_vec,G_d,R_s,theta_s,v_s); %Bias is error squared
Doppler_Bias_Store_SS_8(loop_index,SNR_index)=Doppler_Bias;

%%% Plotting %%%
% Est_Surf_loc=buffer_loc;
% figure(7)
% imagesc(x_vec,y_vec,Est_Surf_loc)
% xlabel('x');
% ylabel('y');

% Est_Surf_Doppler=est_Doppler_signal;
% figure(8)
% plot(Doppler_vec,Est_Surf_Doppler)
% xlabel('Doppler');
% ylabel('Signal Level');
end
end

%% Calculating RMSE
for s=1:SNR_length
MSE_loc(1,s)        = sum(loc_Bias_Store_CS_1(:,s))/(Num_Targets*Loop_Num);
MSE_loc(2,s)        = sum(loc_Bias_Store_CS_4(:,s))/(Num_Targets*Loop_Num);
MSE_loc(3,s)        = sum(loc_Bias_Store_CS_8(:,s))/(Num_Targets*Loop_Num);
MSE_loc(4,s)        = sum(loc_Bias_Store_SS_8(:,s))/(Num_Targets*Loop_Num);

MSE_Doppler(1,s)    = sum(Doppler_Bias_Store_CS_1(:,s))/(Num_Targets*Loop_Num); 
MSE_Doppler(2,s)    = sum(Doppler_Bias_Store_CS_4(:,s))/(Num_Targets*Loop_Num); 
MSE_Doppler(3,s)    = sum(Doppler_Bias_Store_CS_8(:,s))/(Num_Targets*Loop_Num); 
MSE_Doppler(4,s)    = sum(Doppler_Bias_Store_SS_8(:,s))/(Num_Targets*Loop_Num); 
end

RMSE_loc       = sqrt(MSE_loc);
RMSE_Doppler   = sqrt(MSE_Doppler);


%% Saving the data to a file
%save Monte_Carlo_Results.mat
toc
