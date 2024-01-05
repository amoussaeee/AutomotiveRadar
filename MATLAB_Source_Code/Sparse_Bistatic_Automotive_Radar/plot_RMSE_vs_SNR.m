% This script is printing high quality figures.
% Author: Ali Moussa - Last modified: 05/01/2024

%% Loading the data
load("Monte_Carlo_Results.mat")
RMSE_location_Q_1=RMSE_loc(1,:);
RMSE_location_Q_4=RMSE_loc(2,:);
RMSE_location_Q_8=RMSE_loc(3,:);
RMSE_location_MUSIC=RMSE_loc(4,:);

RMSE_Doppler_Q_1=RMSE_Doppler(1,:);
RMSE_Doppler_Q_4=RMSE_Doppler(2,:);
RMSE_Doppler_Q_8=RMSE_Doppler(3,:);
RMSE_Doppler_MUSIC=RMSE_Doppler(4,:);


%% Plotting
figure()
subplot(2,1,1)
plot(SNR_vector,RMSE_location_Q_1,'--x',SNR_vector,RMSE_location_Q_4,'--o',SNR_vector,RMSE_location_Q_8,'--*',SNR_vector,RMSE_location_MUSIC,'--+','MarkerSize',10)


set ( gca, 'YScale', 'log')

legend('Sparse with $Q=1$','Sparse with $Q=4$','Sparse with $Q=8$','MUSIC with $Q=8$','Interpreter','latex') 

xlabel('SNR [dB]','Interpreter','latex')
ylabel('RMSE of location [m]','Interpreter','latex')
xlim([128 145])
grid on
set(gca,'GridLineStyle','--')
set(gca,'FontSize',18) 
orient tall
print -deps RMSE_loc.eps
 
%%

figure()
subplot(2,1,1)
plot(SNR_vector,RMSE_Doppler_Q_1,'--x',SNR_vector,RMSE_Doppler_Q_4,'--o',SNR_vector,RMSE_Doppler_Q_8,'--*',SNR_vector,RMSE_Doppler_MUSIC,'--+','MarkerSize',10)


set ( gca, 'YScale', 'log')

legend('Sparse with $Q=1$','Sparse with $Q=4$','Sparse with $Q=8$','MUSIC with $Q=8$','Interpreter','latex') 

xlabel('SNR [dB]','Interpreter','latex')
ylabel('RMSE of Doppler [m/s]','Interpreter','latex')
xlim([128 145])
grid on
set(gca,'GridLineStyle','--')
set(gca,'FontSize',18) 
orient tall
print -deps RMSE_Doppler.eps


