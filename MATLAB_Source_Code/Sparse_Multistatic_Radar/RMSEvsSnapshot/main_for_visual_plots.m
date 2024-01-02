%% Introduction
% This program is a simulation of an automotive multistatic scenario where 
% the 2D cartesian location is estimated, then Doppler is estimated before
% finally pairing the location and Doppler parameters for K targets. The
% method proposed employs the concept of group sparsity using l21
% minimisation which is solved using the CVX package.

%% Initialisation
clear all;
close all
format long;
format compact;
tic
fprintf('Program executed at:')
fprintf(string(datetime))
%rng(2,'philox')       
%% Define Chirp Parameters
fc=77e09;       %centre frequency
B=150e6;        %chirp bandwidth 
Tc=30e-6;       %modulation period 
c=3e08;         %speed of light
mu=B/Tc;        %modulation rate
K=4;            %number of targets

%% Define  Clock Settings
fs=5e6;         %sampling frequency (~bandwidth in complex sampling)
Ts=1/fs;        %sampling time
N=floor(Tc/Ts); %number of samples in each pulse

%% Maximum Unambiguous Range
Rm=c*fs/mu;     %proposed max detectable range

%% Define ULA Parameters
lambda=c/fc;    %wavelength
d=lambda/2;       %antenna spacing
L=8;            %number of antennas


%% Define Group Sparsity Paramters
Q_vector    = [ 2 4 6 8 ];
Q_length=length(Q_vector);

%% Define Slow Time Parameters
M=128;          %number of chirps   
T=35e-6;        %pulse repetition period
v_min=25; %m/s
v_max=35; %m/s
%% Searching Vehicle and Roadside Sensor Init
%Sensor Init
R_h1=30;
theta_h1=-7.662;
R_h2=30.32941;
theta_h2=11.41;
v_s=25;%m/s


% 2D rectangular grid 
I=21;
J=21;
x_vec=linspace(-4,6,J);
y_vec=linspace(55,65,I);
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
    pol_grid_theta(g)=asind(buffer_x(g)/sqrt(buffer_x(g)^2+buffer_y(g)^2));%Automatically 0 deg when x=0 m.
end
pol_grid=[pol_grid_R',pol_grid_theta'];



%% Define received signal power and noise power
EIRP=1.995262315; %Watts (eq to 33dBm)
Gr=10^(16/10);%linear
RCS=1;

%% Signal Model
q_1k=@(u,a)sqrt((EIRP*Gr*(lambda^2)*RCS)/((4*pi)^3*(u^2)*(sqrt(u^2+R_h1^2-2*u*R_h1*cosd(abs(theta_h1-a)))^2)))...
    *exp(-1i*2*pi*fc*(sqrt(u^2+R_h1^2-2*u*R_h1*cosd(abs(theta_h1-a)))+u)/c);
p_1k=@(u,a)kron(exp((-1i*2*pi*mu*(sqrt(u^2+R_h1^2-2*u*R_h1*cosd(abs(theta_h1-a)))+u)*(0:N-1))/(c*fs)).',exp(-1i*2*pi*fc*d*sind(a)*(0:L-1)/c).'); 

q_2k=@(u,a)sqrt((EIRP*Gr*(lambda^2)*RCS)/((4*pi)^3*(u^2)*(sqrt(u^2+R_h2^2-2*u*R_h2*cosd(abs(theta_h2-a)))^2)))...
    *exp(-1i*2*pi*fc*(sqrt(u^2+R_h2^2-2*u*R_h2*cosd(abs(theta_h2-a)))+u)/c);
p_2k=@(u,a)kron(exp((-1i*2*pi*mu*(sqrt(u^2+R_h2^2-2*u*R_h2*cosd(abs(theta_h2-a)))+u)*(0:N-1))/(c*fs)).',exp(-1i*2*pi*fc*d*sind(a)*(0:L-1)/c).'); 

v_k=@(v) exp(-1i*2*pi*fc*v*(0:M-1)*T/c);

%% Constructing a LNxG location steering matrix
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

%% Constructing a MxG_d Doppler steering matrix
G_d=256;
V_min=min(V_min_var1);
V_max=max(V_max_var1);
Doppler_vec=linspace(V_min,V_max,G_d);
velocity_vec=linspace(v_min,v_max,G_d);
step_Doppler=abs(Doppler_vec(2)-Doppler_vec(1));
V_g=zeros(M,G_d);

for g=1:G_d
   V_g(:,g)=v_k(Doppler_vec(g));
end


%% SNR
SNR=150; %Transmitted
eb_vector_loc        = [0.000015 0.0000300 0.0000500 0.000077]; %0.000077 is good
eb_vector_Doppler_1  = [0.000003 0.0000037 0.0000043 0.0000035];
eb_vector_Doppler_2  = [0.000003 0.0000037 0.0000043 0.0000050];


solver_record_loc=[];
solver_record_Doppler=[];



%% Searching Vehicles
% Target Init
R_1k_search=sqrt(pol_grid_R.^2+R_h1^2-2*pol_grid_R*R_h1.*cosd(theta_h1-pol_grid_theta))+pol_grid_R;
buffer_R_search=sort(R_1k_search);
R_1k_vec=[buffer_R_search(88:88:352)];

buffer_theta=sort(pol_grid_theta);
theta_k_vec=[buffer_theta(88:88:352)];

R_k_vec=(R_1k_vec.^2-R_h1^2)./(2*R_1k_vec-2*R_h1.*cosd(theta_h1-theta_k_vec));
[x_k_vec,y_k_vec]=pol2car(R_k_vec,theta_k_vec);                         


V_1k_vec=[Doppler_vec(26) Doppler_vec(100) Doppler_vec(140) Doppler_vec(200)];    %Doppler of tg from sensor 1
v_k_vec=Doppler2velocity(R_k_vec,theta_k_vec,R_h1,theta_h1,V_1k_vec,v_s);       %velocity of tg 1
V_2k_vec=velocity2Doppler(R_k_vec,theta_k_vec,R_h2,theta_h2,v_k_vec,v_s);       %Doppler of tg  from sensor 2
w_1k_vec=d*sind(theta_k_vec)-T*V_1k_vec;
w_2k_vec=d*sind(theta_k_vec)-T*V_2k_vec;

Target_Positions_pol=[R_k_vec' theta_k_vec'];
Target_Positions_car=[x_k_vec' y_k_vec'];
Target_Doppler=[V_1k_vec' V_2k_vec'];
Target_rav=[Target_Positions_pol,v_k_vec'];
Num_Targets = K;


%% Generate a Noisy Raw Signal using the Signal Model and Target Init
Z_1_n=zeros(L*N,M);
Z_2_n=zeros(L*N,M);

% Signal
for k=1:K
Z_1_n=Z_1_n+q_1k(R_k_vec(k),theta_k_vec(k))*p_1k(R_k_vec(k),theta_k_vec(k))*v_k(V_1k_vec(k));
Z_2_n=Z_2_n+q_2k(R_k_vec(k),theta_k_vec(k))*p_2k(R_k_vec(k),theta_k_vec(k))*v_k(V_2k_vec(k));
end

A_Noise = sqrt(EIRP^2/(10^(SNR/10))/2);
% Noise

W_1 = A_Noise*randn(L*N, M) + A_Noise*1i*randn(L*N, M);
W_2 = A_Noise*randn(L*N, M) + A_Noise*1i*randn(L*N, M);

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

SNR_test=snr((Z_1_n),(W_1)); %Received SNR

%% Snapshot Loop
%for Q_index = 1:Q_length
for Q_index = 4

    Q = Q_vector(Q_index);
    eb_loc      = eb_vector_loc(:,Q_index);
    eb_Doppler_1  = eb_vector_Doppler_1(:,Q_index);
    eb_Doppler_2  = eb_vector_Doppler_2(:,Q_index);
    

%% Test for adjusting the error bounds

noise_test_loc=norm([real(W_1(:,1:Q)),real(W_2(:,1:Q))],'fro');
noise_test_Doppler_1=norm(real(W_1_t(:,1:Q*2)),'fro');
noise_test_Doppler_2=norm(real(W_2_t(:,1:Q*2)),'fro');


%% Proposed %%
% Solve the location  problem using l21 minimisation with M_p=8;
[est_location_signal,~]=l21_minimisation_multistatic([Z_1_comp(:,1:Q),Z_2_comp(:,1:Q)],P_g1,P_g2,eb_loc,G_p,Q);
buffer_loc=reshape(est_location_signal,I,J);
Est_Surf_loc=abs(buffer_loc);

figure(1)
subplot(2,1,1)
imagesc(x_vec,y_vec,Est_Surf_loc/max(Est_Surf_loc(:)))
hold on
h=scatter(x_k_vec(1),y_k_vec(1),100,'red'); 
scatter(x_k_vec(2),y_k_vec(2),100,'red'); 
scatter(x_k_vec(3),y_k_vec(3),100,'red'); 
scatter(x_k_vec(4),y_k_vec(4),100,'red'); 
legend(h,'True Value','Interpreter','latex')
handle=colorbar;
ylabel(handle,'Normalised spectrum level','Interpreter','latex');
set(gca,'FontSize',18) 
set(handle,'TickLabelInterpreter','latex');
xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
grid on
set(gca,'GridLineStyle','--')
set(gca,'TickLabelInterpreter','latex')
hold off
orient tall
print -deps -r300 loc_spec_CS.eps

%% MUSIC %%
% Sensor 1
est_location_signal_1=MUSIC(Z_1_comp(:,1:Q),P_g1,Num_Targets,G_p);
buffer_loc_1=reshape(est_location_signal_1,I,J);


Est_Surf_loc_1=abs(buffer_loc_1);
figure(3)
subplot(2,1,1)
imagesc(x_vec,y_vec,Est_Surf_loc_1/max(Est_Surf_loc_1(:)))
hold on
h=scatter(x_k_vec(1),y_k_vec(1),100,'red'); 
scatter(x_k_vec(2),y_k_vec(2),100,'red'); 
scatter(x_k_vec(3),y_k_vec(3),100,'red'); 
scatter(x_k_vec(4),y_k_vec(4),100,'red'); 
legend(h,'True Value','Interpreter','latex')
handle=colorbar;
ylabel(handle,'Normalised spectrum level','Interpreter','latex');
set(gca,'FontSize',18) 
set(handle,'TickLabelInterpreter','latex');

xlabel('x [m]','Interpreter','latex');
ylabel('y [m]','Interpreter','latex');
grid on
set(gca,'GridLineStyle','--')
set(gca,'TickLabelInterpreter','latex')
hold off
orient tall
print -deps -r300 loc_spec_MUSIC.eps
% 




%% Doppler
est_Doppler_signal_SS1=MUSIC(Z_1_t_comp(:,1:Q*2),V_g,Num_Targets,G_d);
Est_Surf_Doppler_1=abs(est_Doppler_signal_SS1);
[est_Doppler_signal_CS1,~]=l21_minimisation(Z_1_t_comp(:,1:Q*2),V_g,eb_Doppler_1,G_d,2*Q);
Est_Surf_Doppler_CS1=abs(est_Doppler_signal_CS1);


figure(4)
subplot(2,1,1)
holder1=plot(Doppler_vec,Est_Surf_Doppler_CS1/max(Est_Surf_Doppler_CS1(:)),'b-','LineWidth',1);
hold on
holder2=plot(Doppler_vec,Est_Surf_Doppler_1/max(Est_Surf_Doppler_1(:)),'r--','LineWidth',1);
holder3=xline(V_1k_vec,':k');
legend([holder3(1),holder1,holder2],'True Value','GS','MUSIC','Interpreter','latex')
xlabel('Bistatic velocity [m/s]','Interpreter','latex');
ylabel('Normalised spectrum level','Interpreter','latex');
grid on
set(gca,'GridLineStyle','--')
set(gca,'TickLabelInterpreter','latex')
hold off
orient tall
print -deps -r300 Doppler_spec.eps

% % Sensor 2
% est_location_signal_SS2=MUSIC(Z_2_comp(:,1:Q),P_g2,Num_Targets,G_p);
% buffer_loc_2=reshape(est_location_signal_SS2,I,J);
% 
% est_Doppler_signal_2=MUSIC(Z_2_t_comp(:,1:Q*2),V_g,Num_Targets,G_d);
% 


% Est_Surf_loc_2=abs(buffer_loc_2);
% figure(2)
% imagesc(x_vec,y_vec,Est_Surf_loc_2)
% xlabel('x');
% ylabel('y');
% Est_Surf_Doppler_2=abs(est_Doppler_signal_2);
% figure(4)
% plot(Doppler_vec,Est_Surf_Doppler_2)
% xlabel('Doppler');
% ylabel('Signal Level');


end

toc