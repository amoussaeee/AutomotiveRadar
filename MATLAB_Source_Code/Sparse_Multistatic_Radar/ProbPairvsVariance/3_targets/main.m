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
tic
fprintf('Program executed at:')
fprintf(string(datetime))
 %rng(32,'philox')       
       
%%% Define Chirp Parameters %%
fc=77e09;       %centre frequency
B=150e6;        %chirp bandwidth 
Tc=30e-6;       %modulation period 
c=3e08;         %speed of light
mu=B/Tc;        %modulation rate
K=3;            %number of targets

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



%%% Define Slow Time Parameters
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
    pol_grid_theta(g)=asind(buffer_x(g)/sqrt(buffer_x(g)^2+buffer_y(g)^2));%Automatically 0 deg when x=0 m.
end
pol_grid=[pol_grid_R',pol_grid_theta'];



%%% Define received signal power and noise power %%%
EIRP=1.995262315; %Watts (eq to 33dBm)
Gr=10^(16/10);%linear
RCS=1;

%% Signal Model %%
q_1k=@(u,a)sqrt((EIRP*Gr*(lambda^2)*RCS)/((4*pi)^3*(u^2)*(sqrt(u^2+R_h1^2-2*u*R_h1*cosd(abs(theta_h1-a)))^2)))...
    *exp(-1i*2*pi*fc*(sqrt(u^2+R_h1^2-2*u*R_h1*cosd(abs(theta_h1-a)))+u)/c);
p_1k=@(u,a)kron(exp((-1i*2*pi*mu*(sqrt(u^2+R_h1^2-2*u*R_h1*cosd(abs(theta_h1-a)))+u)*(0:N-1))/(c*fs)).',exp(-1i*2*pi*fc*d*sind(a)*(0:L-1)/c).'); 
v_k=@(v) exp(-1i*2*pi*fc*v*(0:M-1)*T/c);

q_2k=@(u,a)sqrt((EIRP*Gr*(lambda^2)*RCS)/((4*pi)^3*(u^2)*(sqrt(u^2+R_h2^2-2*u*R_h2*cosd(abs(theta_h2-a)))^2)))...
    *exp(-1i*2*pi*fc*(sqrt(u^2+R_h2^2-2*u*R_h2*cosd(abs(theta_h2-a)))+u)/c);
p_2k=@(u,a)kron(exp((-1i*2*pi*mu*(sqrt(u^2+R_h2^2-2*u*R_h2*cosd(abs(theta_h2-a)))+u)*(0:N-1))/(c*fs)).',exp(-1i*2*pi*fc*d*sind(a)*(0:L-1)/c).'); 

%%% Constructing a LNxG location steering matrix %%%
P_g1=zeros(L*N,G_p);
V_min_var1=[];
V_max_var1=[];
P_g2=zeros(L*N,G_p);
V_min_var2=[];
V_max_var2=[];
for g=1:G_p
  P_g1(:,g)=p_1k(pol_grid_R(g),pol_grid_theta(g));
  P_g2(:,g)=p_2k(pol_grid_R(g),pol_grid_theta(g));
  V_min_var1=[V_min_var1;velocity2Doppler(pol_grid_R(g),pol_grid_theta(g),R_h1,theta_h1,v_min,v_s)];
  V_max_var1=[V_max_var1;velocity2Doppler(pol_grid_R(g),pol_grid_theta(g),R_h1,theta_h1,v_max,v_s)];
 % V_min_var2=[V_min_var2;velocity2Doppler(pol_grid_R(g),pol_grid_theta(g),R_h2,theta_h2,v_min,v_s)];
 % V_max_var2=[V_max_var2;velocity2Doppler(pol_grid_R(g),pol_grid_theta(g),R_h2,theta_h2,v_max,v_s)];
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
SNR_length=length(SNR_vector);

%%% SD %%%
SD_vector= 0:0.1:0.5;
SD_length= length(SD_vector);

%%% Loops %%%
Loop_Num=1000;

prob_correct_pairing_LS_1=zeros(SD_length,SNR_length);
prob_correct_pairing_CC_1=zeros(SD_length,SNR_length);

pairing_record_LS_1=strings(Loop_Num,SNR_length);
pairing_record_CC_1=strings(Loop_Num,SNR_length);




%% Searching Vehicle
% Target Init
x_k_vec=[x_vec(3) x_vec(7) x_vec(11)];                          %X of target 
y_k_vec=[y_vec(3) y_vec(7) y_vec(11)];                        %Y of target 
[R_k_vec,theta_k_vec] = car2pol(x_k_vec,y_k_vec);

V_1k_vec=[Doppler_vec(20)  Doppler_vec(60) Doppler_vec(110)];    %Doppler of tg from sensor 1
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
%% Standard Deviation Loop Start %%
for SD_index=1:SD_length
%for SD_index=1

    SD=SD_vector(SD_index);

   
%% Monte Carlo Loop Start %%
for loop_index = 1:Loop_Num
    fprintf('\n Current Monte Carlo Trial Number:')
    fprintf(string(loop_index))

Target_rav_est=normrnd(Target_rav,SD);
est_R_k       = Target_rav_est(:,1);
est_theta_k   = Target_rav_est(:,2);
est_Doppler_k = velocity2Doppler(Target_rav_est(:,1),Target_rav_est(:,2),R_h1,theta_h1,Target_rav_est(:,3),v_s);
W_1_n = randn(L*N, M) + 1i*randn(L*N, M);
W_2_n = randn(L*N, M) + 1i*randn(L*N, M);

%% SNR Loop Start %%
 %for SNR_index = 1:SNR_length
 for SNR_index = 5

    SNR=SNR_vector(SNR_index);
A_Noise = sqrt(EIRP^2/(10^(SNR/10))/2);

W_1 = A_Noise*W_1_n;
W_2 = A_Noise*W_2_n;

Z_1 = Z_1_n+W_1;
Z_2 = Z_2_n+W_2;
Z_1_t=Z_1.';
Z_2_t=Z_2.';


%% SNR test
SNR_test=snr((Z_1_n),(W_1));

%% Parameter pairing %%
%%% Cross Correlation Method %%%
cov_length=L;
N_vec=((1:N)-1)*L+1;
Z_a=Z_1(:,1);
Z_a=reshape(Z_a,[L N]);
Z_a=Z_a(1:cov_length,:);
Z_v=Z_1_t(1:cov_length,N_vec);
z_av=zeros(cov_length,1);
for n=1:N
z_av=z_av+(Z_a(:,n).*conj(Z_v(:,n)));
end
z_av=z_av/N;
Rcc=toeplitz(z_av,conj(z_av));
W_vec=linspace(-3.85e-4,-1.3039e-3,256);
av_k=@(v) exp(-1i*2*pi*fc*v*(0:cov_length-1)/c).';
av_g=zeros(cov_length,length(W_vec));
for g=1:length(W_vec)
av_g(:,g)=av_k(W_vec(g));
end

[W_k,~,~]=esprit_propagator_music(Rcc,av_g,cov_length,Num_Targets,fc,length(W_vec));
Target_comb=[];
perm=[];
for i=1:K
    for j=1:K
        perm=[perm; [i j]];
        Target_comb=[Target_comb;[est_theta_k(i) est_Doppler_k(j)]];
    end
end
pairing_fun=zeros(K*K,K);
min_comb=zeros(K,1);
for k=1:K
for comb=1:K*K
pairing_fun(comb,k)=norm(exp(-1i*2*pi*fc*W_k(k)/c)-exp(-1i*2*pi*fc*(d*sind(Target_comb(comb,1))-T*Target_comb(comb,2))/c));
end
[~,min_comb(k)]=min(pairing_fun(:,k),[],1);
end
pairing_test=ismember(1,min_comb) & ismember(5,min_comb) & ismember(9,min_comb);
if pairing_test==1
    pairing_record_CC_1(loop_index,SNR_index)="success";
else 
    pairing_record_CC_1(loop_index,SNR_index)="fail";
end

%%% Least Squares Method %%%
p_k_est=zeros(L*N,K);
V_1k_est=zeros(M,K);
V_2k_est=zeros(M,K);
matched_steering=zeros(L*N*M,K);
perms_LS=perms(1:K);
pair_norm=zeros(factorial(K),1);
% pair_norm_test=zeros(factorial(K),1);
for p=1:factorial(K)
    Doppler=est_Doppler_k(perms_LS(p,:));
    for k=1:K
        matched_steering(:,k)=kron(v_k(Doppler(k)).',p_1k(est_R_k(k),est_theta_k(k)));
    end
   % x_pair=lsqr((matched_steering),(Z_1(:)));
   % pair_norm_test(p)=norm((Z_1(:))-(matched_steering)*x_pair);
   pair_norm(p)=norm((Z_1(:))-(matched_steering)*((matched_steering'*matched_steering)\matched_steering'*Z_1(:)));

end
[~,min_perm]=min(pair_norm,[],1);
if min_perm==factorial(K)
    pairing_record_LS_1(loop_index,SNR_index)="success";
else 
    pairing_record_LS_1(loop_index,SNR_index)="fail";
end

end

end

for s=1:SNR_length
prob_correct_pairing_LS_1(SD_index,s)    = sum(count(pairing_record_LS_1(:,s),"success"))/Loop_Num;
prob_correct_pairing_CC_1(SD_index,s)    = sum(count(pairing_record_CC_1(:,s),"success"))/Loop_Num;
end

end

%% Saving the data to a file %%
save ProbvsNum_test_1000.mat
toc