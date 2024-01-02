%% Introduction
% This program is a simulation of an automotive multistatic scenario where 
% the 2D cartesian location is estimated, then Doppler is estimated before
% finally pairing the location and Doppler parameters for K targets. The
% method proposed employs the concept of group sparsity using l21
% minimisation which is solved using the CVX package.

%% Initialisation %%
clear all;
close all
format long;
format compact;
%tic
fprintf('Program executed at:')
fprintf(string(datetime))
%rng(32,'philox')  


%%% Define Chirp Parameters %%
fc=77e09;       %centre frequency
B=150e6;        %chirp bandwidth 
Tc=30e-6;       %modulation period 
c=3e08;         %speed of light
mu=B/Tc;        %modulation rate
K=1;            %number of targets

%%% Define  Clock Settings %%%
fs=5e6;         %sampling frequency (~bandwidth in complex sampling)
Ts=1/fs;        %sampling time
N=floor(Tc/Ts); %number of samples in each pulse

%%% Maximum Unambiguous Range %%%
Rm=c*fs/mu;     %proposed max detectable range

%%% Define ULA Parameters %%%
lambda=c/fc;    %wavelength
d=lambda/2;       %antenna spacing
L=8;            %number of antennas


%%% Define Group Sparsity Paramters %%%
% Location Estimation
Q_l = 8;       %number of processed pulses

% Doppler Estimation
Q_d = 16;      %number of processed snapshots

%%% Define Slow Time Parameters %%%
M=128;          %number of chirps   
T=35e-6;        %pulse repetition period
v_min=25; %m/s
v_max=35; %m/s

%%% Searching Vehicle and Roadside Sensor Init %%%
% Sensor Init
R_h1=30;
theta_h1=-7.662;
R_h2=30.32941;
theta_h2=11.41;
v_s=25;%m/s


%%% 2D rectangular grid %%%
I=13;
J=13;
x_vec=linspace(0,6,J);
y_vec=linspace(57,63,I);
step_x=abs(x_vec(2)-x_vec(1));
step_y=abs(y_vec(2)-y_vec(1));

rect_grid_x=zeros(I,J);
rect_grid_y=zeros(I,J);
for j=1:J
    rect_grid_x(:,j)=x_vec(j);
end

for i=1:I
    rect_grid_y(i,:)=y_vec(i);
end

% Convert 2D rectangular grid to 1D polar grid of length G=IxJ
G_p=I*J;
pol_grid_R=zeros(1,G_p);
pol_grid_theta=zeros(1,G_p);
buffer_x=rect_grid_x(:);
buffer_y=rect_grid_y(:);

for g=1:G_p
    pol_grid_R(g)=sqrt(buffer_x(g)^2+buffer_y(g)^2);
    pol_grid_theta(g)=asind(buffer_x(g)/sqrt(buffer_x(g)^2+buffer_y(g)^2)); % Automatically 0 deg when x=0 m.
end
pol_grid=[pol_grid_R',pol_grid_theta'];

loc=1:K;
dop=1:K;

%%% Define received signal power and noise power %%%
EIRP=1.995262315; %Watts (eq to 33dBm)
Gr=10^(16/10);%linear
RCS=1;

%% Signal Model %%
q_1k=@(u,a)sqrt((EIRP*Gr*(lambda^2)*RCS)/((4*pi)^3*(u^2)*(sqrt(u^2+R_h1^2-2*u*R_h1*cosd(abs(theta_h1-a)))^2)))...
    *exp(-1i*2*pi*fc*(sqrt(u^2+R_h1^2-2*u*R_h1*cosd(abs(theta_h1-a)))+u)/c);
q_2k=@(u,a)sqrt((EIRP*Gr*(lambda^2)*RCS)/((4*pi)^3*(u^2)*(sqrt(u^2+R_h2^2-2*u*R_h2*cosd(abs(theta_h2-a)))^2)))...
    *exp(-1i*2*pi*fc*(sqrt(u^2+R_h2^2-2*u*R_h2*cosd(abs(theta_h2-a)))+u)/c);

p_1k=@(u,a)kron(exp((-1i*2*pi*mu*(sqrt(u^2+R_h1^2-2*u*R_h1*cosd(abs(theta_h1-a)))+u)*(0:N-1))/(c*fs)).',exp(-1i*2*pi*fc*d*sind(a)*(0:L-1)/c).'); 
p_2k=@(u,a)kron(exp((-1i*2*pi*mu*(sqrt(u^2+R_h2^2-2*u*R_h2*cosd(abs(theta_h2-a)))+u)*(0:N-1))/(c*fs)).',exp(-1i*2*pi*fc*d*sind(a)*(0:L-1)/c).'); 

v_k=@(v) exp(-1i*2*pi*fc*v*(0:M-1)*T/c);



%%% Constructing a LNxG location steering matrix %%%
P_g1=zeros(L*N,G_p);
P_g2=zeros(L*N,G_p);

V_min_var1=[];
V_max_var1=[];

for g=1:G_p
  P_g1(:,g)=p_1k(pol_grid_R(g),pol_grid_theta(g));
  P_g2(:,g)=p_2k(pol_grid_R(g),pol_grid_theta(g));
  V_min_var1=[V_min_var1;velocity2Doppler(pol_grid_R(g),pol_grid_theta(g),R_h1,theta_h1,v_min,v_s)];
  V_max_var1=[V_max_var1;velocity2Doppler(pol_grid_R(g),pol_grid_theta(g),R_h1,theta_h1,v_max,v_s)];
end

%%% Constructing a MxG_d Doppler steering matrix %%%
G_d=128;
V_min=min(V_min_var1);
V_max=max(V_max_var1);
Doppler_vec=linspace(V_min,V_max,G_d);
velocity_vec=linspace(v_min,v_max,G_d);
step_Doppler=abs(Doppler_vec(2)-Doppler_vec(1));
V_g=zeros(M,G_d);
for g=1:G_d
   V_g(:,g)=v_k(Doppler_vec(g));
end


%%% SNR %%%
SNR_vector=120:10:160;

eb_vector_loc        = [0.000198 0.000070 0.000040 0.000035 0.000034];
eb_vector_Doppler_1  = [0.000065 0.000020 0.000007 0.000002 0.000001];
eb_vector_Doppler_2  = [0.000065 0.000020 0.000007 0.000002 0.000001];


SNR_length=length(SNR_vector);

%%% Loops %%%
Loop_Num=1;
Methods_num=4;
MSE_loc=zeros(Methods_num,SNR_length);
MSE_velocity=zeros(Methods_num,SNR_length);
MSE_Doppler_1=zeros(2,SNR_length);
MSE_Doppler_2=zeros(2,SNR_length);


% Compressive Sensing buffers
loc_Bias_Store_CS         = zeros(K*Loop_Num,SNR_length);
velocity_Bias_Store_CS    = zeros(K*Loop_Num,SNR_length);
Doppler_Bias_Store_CS_1   = zeros(K*Loop_Num,SNR_length);
Doppler_Bias_Store_CS_2   = zeros(K*Loop_Num,SNR_length);



solver_record_loc=[];
solver_record_Doppler=[];

% MUSIC buffers
loc_Bias_Store_SS         = zeros(K*Loop_Num,SNR_length);
velocity_Bias_Store_SS    = zeros(K*Loop_Num,SNR_length);
Doppler_Bias_Store_SS_1   = zeros(K*Loop_Num,SNR_length);
Doppler_Bias_Store_SS_2   = zeros(K*Loop_Num,SNR_length);


%% Monte Carlo Loop Start %%
for loop_index = 1:Loop_Num
    fprintf('Current Monte Carlo Trial Number:')
    fprintf(string(loop_index))
%% Searching Vehicle
% Target Init
x_k1=unifrnd(x_vec(6)+3*step_x/8,x_vec(7)-3*step_x/8);                            % X of tg 1
y_k1=unifrnd(y_vec(6)+3*step_y/8,y_vec(7)-3*step_y/8);                            % Y of tg 1

%x_k1=(x_vec(6)+x_vec(7))/2;                          % X of tg 1
%y_k1=(y_vec(6)+y_vec(7))/2;                          % Y of t 1


V_1k1=unifrnd(Doppler_vec(72)+3*step_Doppler/8,Doppler_vec(73)-3*step_Doppler/8); % Doppler of tg 1 from sensor 1

%V_1k1=(Doppler_vec(72)+Doppler_vec(73))/2; % Doppler of tg 1 from sensor 1


[R_k1,theta_k1]= car2pol(x_k1,y_k1);                                          % range and DOA of tg 1
v_k1=Doppler2velocity(R_k1,theta_k1,R_h1,theta_h1,V_1k1,v_s);                 % velocity of tg 1
V_2k1=velocity2Doppler(R_k1,theta_k1,R_h2,theta_h2,v_k1,v_s);

% Bistatic angle based on target 1 and sensor 1
%beta_1k1=acosd((R_k1^2+R_h1^2-2*R_k1*R_h1*cosd(abs(theta_h1-theta_k1))+R_k1^2-R_h1^2)/(2*sqrt(R_k1^2+R_h1^2-2*R_k1*R_h1*cosd(abs(theta_h1-theta_k1)))*R_k1));

% Radial range resolution based on target 1 and sensor 1
%range_res=c/(2*B*cosd(beta_1k1/2));

% DOA resolution based on target 1
%DOA_res=c/(fc*L*d*cosd(theta_k1));

% Cartesian resolution
%x_res=range_res*sind(DOA_res);
%y_res=range_res*cosd(DOA_res);



Target_Positions_pol=[R_k1 theta_k1];
Target_Positions_car=[x_k1 y_k1];
Target_Doppler=[V_1k1 V_2k1];
Target_rav=[Target_Positions_pol,v_k1];
Num_Targets = K;

%% Generate a Noisy Raw Signal using the Signal Model and Target Init

% Signal
Z_1_n=q_1k(R_k1,theta_k1)*p_1k(R_k1,theta_k1)*v_k(V_1k1);
Z_2_n=q_2k(R_k1,theta_k1)*p_2k(R_k1,theta_k1)*v_k(V_2k1);


    
W_1_n = randn(L*N, M) + 1i*randn(L*N, M);
W_2_n = randn(L*N, M) + 1i*randn(L*N, M);

%% SNR Loop %%
%for SNR_index = 1:SNR_length
for SNR_index = 5

    SNR = SNR_vector(SNR_index);
    eb_loc      = eb_vector_loc(:,SNR_index);
    eb_Doppler_1  = eb_vector_Doppler_1(:,SNR_index);
    eb_Doppler_2  = eb_vector_Doppler_2(:,SNR_index);
    A_Noise = sqrt(EIRP^2/(10^(SNR/10))/2);


W_1 = A_Noise*W_1_n;
W_2 = A_Noise*W_2_n;
W_1_t=W_1.';
W_2_t=W_2.';

Z_1_comp=Z_1_n+W_1;
Z_2_comp=Z_2_n+W_2;
Z_1_t_comp=Z_1_comp.';
Z_2_t_comp=Z_2_comp.';

Z_1=real(Z_1_comp);
Z_2=real(Z_2_comp);
Z_1_t=Z_1.';
Z_2_t=Z_2.';

%% Test for adjusting the error bounds
SNR_test=snr((Z_1_n),(W_1));

noise_test_loc       = norm([real(W_1(:,1:Q_l)),real(W_2(:,1:Q_l))],'fro');
noise_test_Doppler_1 = norm(real(W_1_t(:,1:Q_d)),'fro');
noise_test_Doppler_2 = norm(real(W_2_t(:,1:Q_d)),'fro');


%% Proposed %%
%%% Solve the location  problem using l21 minimisation with Q_l=8 %%%
[est_location_signal,solver_status_loc]=l21_minimisation_multistatic([Z_1(:,1:Q_l),Z_2(:,1:Q_l)],P_g1,P_g2,eb_loc,G_p,Q_l);
buffer_loc=reshape(est_location_signal,I,J);

if strcmp(solver_status_loc,'Solved')==0
  solver_record_loc=[solver_record_loc; [solver_status_loc string(SNR) string(eb_loc) string(noise_test_loc) string(loop_index) 'location']];
end

[est_loc_Value, loc_Bias,~] = Location_Estimation_Bias(Num_Targets, buffer_loc, Target_Positions_pol, Target_Positions_car, pol_grid, x_vec, y_vec,I,J); %Bias is error squared
loc_Bias_Store_CS(loop_index,SNR_index)=loc_Bias;
[est_R_k, est_theta_k]=car2pol(est_loc_Value(:,1),est_loc_Value(:,2));

% Est_Surf_loc=buffer_loc;
% figure(1)
% imagesc(x_vec,y_vec,Est_Surf_loc)
% xlabel('x');
% ylabel('y');

%%% Solve the Doppler problem using l21 minimisation with Q_d=16 %%%
% from sensor 1
[est_Doppler_signal_1,solver_status_Doppler]=l21_minimisation(Z_1_t(:,1:Q_d),V_g,eb_Doppler_1,G_d,Q_d);
if strcmp(solver_status_Doppler,'Solved')==0
    solver_record_Doppler=[solver_record_Doppler; [solver_status_Doppler string(SNR) string(eb_Doppler_1) string(noise_test_Doppler_1) string(loop_index) "sensor 1"]];
end
[est_Doppler_Value_1,Doppler_Bias_1,~] = Doppler_Estimation_Bias(Num_Targets, est_Doppler_signal_1, Target_Doppler(:,1), Doppler_vec);
Doppler_Bias_Store_CS_1(loop_index,SNR_index)=Doppler_Bias_1;

% Est_Surf_Doppler_1=est_Doppler_signal_1;
% figure(2)
% plot(Doppler_vec,Est_Surf_Doppler_1)
% xlabel('Doppler');
% ylabel('Signal Level');

%%% Solve the Doppler problem using l21 minimisation with Q_d=16 %%%
% from sensor 2
[est_Doppler_signal_2,solver_status_Doppler]=l21_minimisation(Z_2_t(:,1:Q_d),V_g,eb_Doppler_2,G_d,Q_d);
if strcmp(solver_status_Doppler,'Solved')==0
    solver_record_Doppler=[solver_record_Doppler; [solver_status_Doppler string(SNR) string(eb_Doppler_2) string(noise_test_Doppler_2) string(loop_index) "sensor 2"]];
end
[est_Doppler_Value_2,Doppler_Bias_2,~] = Doppler_Estimation_Bias(Num_Targets, est_Doppler_signal_2, Target_Doppler(:,2), Doppler_vec);
Doppler_Bias_Store_CS_2(loop_index,SNR_index)=Doppler_Bias_2;

% Est_Surf_Doppler_2=est_Doppler_signal_2;
% figure(3)
% plot(Doppler_vec,Est_Surf_Doppler_2)
% xlabel('Doppler');
% ylabel('Signal Level');

%% Calculating velocity and bias from CS %%
est_parameters_1=zeros(K,3);
est_parameters_2=zeros(K,3);

for k=1:K
    est_parameters_1(k,:)=[est_R_k(loc(k)) est_theta_k(loc(k)) est_Doppler_Value_1(dop(k))];
    est_parameters_2(k,:)=[est_R_k(loc(k)) est_theta_k(loc(k)) est_Doppler_Value_2(dop(k))];
end


for k=1:K
    est_parameters_1(k,3)=Doppler2velocity(est_parameters_1(k,1),est_parameters_1(k,2),R_h1,theta_h1,est_parameters_1(k,3),v_s);
    est_parameters_2(k,3)=Doppler2velocity(est_parameters_2(k,1),est_parameters_2(k,2),R_h2,theta_h2,est_parameters_2(k,3),v_s);
end

est_parameters_final=(est_parameters_1+est_parameters_2)/2;
velocity_Bias = (est_parameters_final(:,3)-Target_rav(:,3)).^2;
velocity_Bias_Store_CS(loop_index,SNR_index)=velocity_Bias;


tic
%% MUSIC %%
% Sensor 1
est_location_signal_1=MUSIC(Z_1_comp,(P_g1),Num_Targets,G_p,Q_l);
buffer_loc_1=reshape(est_location_signal_1,I,J);

est_Doppler_signal_1=MUSIC(Z_1_t_comp,(V_g),Num_Targets,G_d,Q_d);
[est_Doppler_Value_1,Doppler_Bias_1,~] = Doppler_Estimation_Bias(Num_Targets, abs(est_Doppler_signal_1), Target_Doppler(:,1), Doppler_vec);
Doppler_Bias_Store_SS_1(loop_index,SNR_index)=Doppler_Bias_1;


% Sensor 2
est_location_signal_2=MUSIC(Z_2_comp,(P_g2),Num_Targets,G_p,Q_l);
buffer_loc_2=reshape(est_location_signal_2,I,J);

est_Doppler_signal_2=MUSIC(Z_2_t_comp,(V_g),Num_Targets,G_d,Q_d);
[est_Doppler_Value_2,Doppler_Bias_2,~] = Doppler_Estimation_Bias(Num_Targets, abs(est_Doppler_signal_2), Target_Doppler(:,2), Doppler_vec);
Doppler_Bias_Store_SS_2(loop_index,SNR_index)=Doppler_Bias_2;


% Location Average
[est_loc_Value, loc_Bias,~] = Location_Estimation_Bias_Average(Num_Targets, abs(buffer_loc_1),abs(buffer_loc_2), Target_Positions_pol, Target_Positions_car, pol_grid, x_vec, y_vec,I,J); 
loc_Bias_Store_SS(loop_index,SNR_index)=loc_Bias;
[est_R_k, est_theta_k]=car2pol(est_loc_Value(:,1),est_loc_Value(:,2));


% Est_Surf_loc_1=abs(buffer_loc_1);
% figure(1)
% imagesc(x_vec,y_vec,Est_Surf_loc_1)
% xlabel('x');
% ylabel('y');
% 
% 
% Est_Surf_loc_2=abs(buffer_loc_2);
% figure(2)
% imagesc(x_vec,y_vec,Est_Surf_loc_2)
% xlabel('x');
% ylabel('y');
% 
% Est_Surf_Doppler_1=abs(est_Doppler_signal_1);
% figure(3)
% plot(Doppler_vec,Est_Surf_Doppler_1)
% xlabel('Doppler');
% ylabel('Signal Level');
% 
% Est_Surf_Doppler_2=abs(est_Doppler_signal_2);
% figure(4)
% plot(Doppler_vec,Est_Surf_Doppler_2)
% xlabel('Doppler');
% ylabel('Signal Level');





%% Calculating velocity and bias from MUSIC %%
est_parameters_1=zeros(K,3);
est_parameters_2=zeros(K,3);

for k=1:K
    est_parameters_1(k,:)=[est_R_k(loc(k)) est_theta_k(loc(k)) est_Doppler_Value_1(dop(k))];
    est_parameters_2(k,:)=[est_R_k(loc(k)) est_theta_k(loc(k)) est_Doppler_Value_2(dop(k))];
end

for k=1:K
    est_parameters_1(k,3)=Doppler2velocity(est_parameters_1(k,1),est_parameters_1(k,2),R_h1,theta_h1,est_parameters_1(k,3),v_s);
    est_parameters_2(k,3)=Doppler2velocity(est_parameters_2(k,1),est_parameters_2(k,2),R_h2,theta_h2,est_parameters_2(k,3),v_s);
end

est_parameters_final=(est_parameters_1+est_parameters_2)/2;
velocity_Bias = (est_parameters_final(:,3)-Target_rav(:,3)).^2;
velocity_Bias_Store_SS(loop_index,SNR_index)=velocity_Bias;
toc
end
end

%% Calculating RMSE %%

  for s=1:SNR_length
MSE_loc(1,s)      = sum(loc_Bias_Store_CS(:,s))/(Num_Targets*Loop_Num);
MSE_loc(2,s)      = sum(loc_Bias_Store_SS(:,s))/(Num_Targets*Loop_Num);
MSE_velocity(1,s) = sum(velocity_Bias_Store_CS(:,s))/(Num_Targets*Loop_Num);
MSE_velocity(2,s) = sum(velocity_Bias_Store_SS(:,s))/(Num_Targets*Loop_Num);


MSE_Doppler_1(1,s) = sum(Doppler_Bias_Store_CS_1(:,s))/(Num_Targets*Loop_Num); 
MSE_Doppler_2(1,s) = sum(Doppler_Bias_Store_CS_2(:,s))/(Num_Targets*Loop_Num); 
MSE_Doppler_1(2,s) = sum(Doppler_Bias_Store_SS_1(:,s))/(Num_Targets*Loop_Num); 
MSE_Doppler_2(2,s) = sum(Doppler_Bias_Store_SS_2(:,s))/(Num_Targets*Loop_Num); 
 
  end

RMSE_loc       = sqrt(MSE_loc);
RMSE_velocity  = sqrt(MSE_velocity);
RMSE_Doppler_1 = sqrt(MSE_Doppler_1);
RMSE_Doppler_2 = sqrt(MSE_Doppler_2);

%% Saving the data to a file %%
%save epsilon_test_30.mat
%toc