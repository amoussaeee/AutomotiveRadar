%% Loading the data
load("RMSEvsSnapshot_test_1000.mat")
RMSE_location_CS=RMSE_loc(1,:);
RMSE_location_MUSIC=RMSE_loc(2,:);

RMSE_Doppler_CS=RMSE_velocity(1,:);
RMSE_Doppler_MUSIC=RMSE_velocity(2,:);



%% Plotting
figure()
subplot(2,1,1)
plot(Q_vector,RMSE_location_CS,'b-o',Q_vector,RMSE_location_MUSIC,'r-x','MarkerSize',10,'LineWidth',1)


set ( gca, 'YScale', 'log')

legend('GS-Joint','MUSIC-Average','Interpreter','latex','Location','northeast') 

xlabel('Number of pulses','Interpreter','latex')
ylabel('RMSE of location [m]','Interpreter','latex')
%xlim([120 160])
grid on
set(gca,'GridLineStyle','--')
set(gca,'FontSize',18) 
set(gca,'TickLabelInterpreter','latex')

orient tall
print -depsc -r300 RMSEvsPulse_loc.eps
 
%%

figure()
subplot(2,1,1)
plot(2*Q_vector,RMSE_Doppler_CS,'b-o',2*Q_vector,RMSE_Doppler_MUSIC,'r-x','MarkerSize',10,'LineWidth',1)


set ( gca, 'YScale', 'log')

legend('GS-Average','MUSIC-Average','Interpreter','latex','Location','northeast') 

xlabel('Number of snapshots','Interpreter','latex')
ylabel('RMSE of velocity [m/s]','Interpreter','latex')
%xlim([120 160])
%yticks([0.11 0.12 0.13 0.14])
%yticklabels({'0.11', '0.12', '0.13', '0.14'})
grid on
set(gca,'GridLineStyle','--')
set(gca,'FontSize',18) 
set(gca,'TickLabelInterpreter','latex')

orient tall
print -depsc -r300 RMSEvsSnapshot_Doppler.eps


