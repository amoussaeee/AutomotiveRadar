%% Loading the data
%load("SNR_test_500.mat")
RMSE_location_CS=RMSE_loc(1,:);
RMSE_location_MUSIC=RMSE_loc(2,:);

RMSE_Doppler_CS=RMSE_velocity(1,:);
RMSE_Doppler_MUSIC=RMSE_velocity(2,:);

SNR_vector=120:10:160;
%% Plotting
figure()
subplot(2,1,1)
plot(SNR_vector,RMSE_location_CS,'--o',SNR_vector,RMSE_location_MUSIC,'--x','MarkerSize',10)


%set ( gca, 'YScale', 'log')

legend('GS-MMV','\emph{Average}-MUSIC-MMV','Interpreter','latex') 

xlabel('SNR [dB]','Interpreter','latex')
ylabel('RMSE of location [m]','Interpreter','latex')
xlim([120 160])
%grid on
%set(gca,'GridLineStyle','--')
set(gca,'FontSize',18) 
orient tall
print -deps RMSEvsSNR_loc.eps
 
%%

figure()
subplot(2,1,1)
plot(SNR_vector,RMSE_Doppler_CS,'--o',SNR_vector,RMSE_Doppler_MUSIC,'--x','MarkerSize',10)


%set ( gca, 'YScale', 'log')

legend('GS-MMV','\emph{Average}-MUSIC-MMV','Interpreter','latex') 

xlabel('SNR [dB]','Interpreter','latex')
ylabel('RMSE of velocity [m/s]','Interpreter','latex')
xlim([120 160])
%grid on
%set(gca,'GridLineStyle','--')
set(gca,'FontSize',18) 
orient tall
print -deps RMSEvsSNR_Doppler.eps


