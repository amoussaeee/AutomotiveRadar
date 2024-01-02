%% Loading the data
load("ProbvsNum_test_2000_100_dB.mat")
Prob_location_LS=prob_correct_pairing_LS_1;
Prob_location_CC=prob_correct_pairing_CC_1;
SNR_vector=100:10:160;
sigma=0:0.1:0.7;
colormap(parula(length(sigma)))
parulacustom = parula(length(sigma));


%% Plotting

figure()
subplot(2,1,1)
lo=plot(SNR_vector,Prob_location_LS(1,:),'k-','MarkerSize',10,'LineWidth',1);hold on
lx=plot(SNR_vector,Prob_location_CC(1,:),'k-.','MarkerSize',10,'LineWidth',1);
L(1)=plot(nan,nan,'ko','MarkerSize',10,'LineWidth',1);
L(2)=plot(nan,nan,'ksquare','MarkerSize',10,'LineWidth',1);
L(3)=plot(nan,nan,'k>','MarkerSize',10,'LineWidth',1);
L(4)=plot(nan,nan,'kx','MarkerSize',10,'LineWidth',1);
L(5)=plot(nan,nan,'k<','MarkerSize',10,'LineWidth',1);
L(6)=plot(nan,nan,'k*','MarkerSize',10,'LineWidth',1);
L(7)=plot(nan,nan,'k+','MarkerSize',10,'LineWidth',1);
L(8)=plot(nan,nan,'k^','MarkerSize',10,'LineWidth',1);

plot(SNR_vector,Prob_location_LS(1,:),'-o','Color',parulacustom(1,:),'MarkerSize',10,'LineWidth',1);
plot(SNR_vector,Prob_location_CC(1,:),'-.o','Color',parulacustom(1,:),'MarkerSize',10,'LineWidth',1);

plot(SNR_vector,Prob_location_LS(2,:),'-square','Color',parulacustom(2,:),'MarkerSize',10,'LineWidth',1);
plot(SNR_vector,Prob_location_CC(2,:),'-.square','Color',parulacustom(2,:),'MarkerSize',10,'LineWidth',1);

plot(SNR_vector,Prob_location_LS(3,:),'->','Color',parulacustom(3,:),'MarkerSize',10,'LineWidth',1);
plot(SNR_vector,Prob_location_CC(3,:),'-.>','Color',parulacustom(3,:),'MarkerSize',10,'LineWidth',1);

plot(SNR_vector,Prob_location_LS(4,:),'-x','Color',parulacustom(4,:),'MarkerSize',10,'LineWidth',1);
plot(SNR_vector,Prob_location_CC(4,:),'-.x','Color',parulacustom(4,:),'MarkerSize',10,'LineWidth',1);

plot(SNR_vector,Prob_location_LS(5,:),'-<','Color',parulacustom(5,:),'MarkerSize',10,'LineWidth',1);
plot(SNR_vector,Prob_location_CC(5,:),'-.<','Color',parulacustom(5,:),'MarkerSize',10,'LineWidth',1);

plot(SNR_vector,Prob_location_LS(6,:),'-*','Color',parulacustom(6,:),'MarkerSize',10,'LineWidth',1);
plot(SNR_vector,Prob_location_CC(6,:),'-.*','Color',parulacustom(6,:),'MarkerSize',10,'LineWidth',1);

plot(SNR_vector,Prob_location_LS(7,:),'-+','Color',parulacustom(7,:),'MarkerSize',10,'LineWidth',1);
plot(SNR_vector,Prob_location_CC(7,:),'-.+','Color',parulacustom(7,:),'MarkerSize',10,'LineWidth',1);

plot(SNR_vector,Prob_location_LS(8,:),'-^','Color',parulacustom(8,:),'MarkerSize',10,'LineWidth',1);
plot(SNR_vector,Prob_location_CC(8,:),'-.^','Color',parulacustom(8,:),'MarkerSize',10,'LineWidth',1);


legend([lo lx L(1) L(2) L(3) L(4) L(5) L(6) L(7) L(8)],'Pair-LS','Pair-CC-ESPRIT','$\sigma_e=0$','$\sigma_e=0.1$','$\sigma_e=0.2$','$\sigma_e=0.3$','$\sigma_e=0.4$','$\sigma_e=0.5$','$\sigma_e=0.6$','$\sigma_e=0.7$','Interpreter','latex','Location','southwest');
xlabel('Input SNR [dB]','Interpreter','latex')
ylabel('Probability of successful pairing','Interpreter','latex')

%xlim([2 6])
% xticklabels([2 3 4 5])
% xticks([2 3 4 5])

grid on
set(gca,'GridLineStyle','--')
set(gca,'FontSize',18) 
set(gca,'TickLabelInterpreter','latex')
%set(gca, 'YScale', 'log')
colormap(parula(length(sigma)))
cb = colorbar;
caxis([0 0.7])
set(cb, 'TickLabelInterpreter', 'latex');
ylabel(cb,'$\sigma_e$','Interpreter','latex','FontSize',18)


orient tall
print -depsc -r300 ProbvsSNR.eps
