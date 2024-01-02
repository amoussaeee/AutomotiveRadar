%% Loading the data
load("ProbvsNum_CC.mat")
load("ProbvsNum_LS.mat")
Prob_location_LS=prob_correct_pairing_LS_1;
Prob_location_CC=prob_correct_pairing_CC_1;
Num_targets=[2 3 4 5 6];
sigma=0:0.1:0.5;
colormap(parula(length(sigma)))
parulacustom = parula(length(sigma));



L=zeros(length(sigma),1);
%% Plotting

figure()
subplot(2,1,1)
lo=plot(Num_targets,Prob_location_LS(1,:),'k-','MarkerSize',10,'LineWidth',1);hold on
lx=plot(Num_targets,Prob_location_CC(1,:),'k-.','MarkerSize',10,'LineWidth',1);
L(1)=plot(nan,nan,'ko','MarkerSize',10,'LineWidth',1);
L(2)=plot(nan,nan,'ksquare','MarkerSize',10,'LineWidth',1);
L(3)=plot(nan,nan,'k>','MarkerSize',10,'LineWidth',1);
L(4)=plot(nan,nan,'kx','MarkerSize',10,'LineWidth',1);
L(5)=plot(nan,nan,'k<','MarkerSize',10,'LineWidth',1);
L(6)=plot(nan,nan,'k*','MarkerSize',10,'LineWidth',1);

plot(Num_targets,Prob_location_LS(1,:),'-o','Color',parulacustom(1,:),'MarkerSize',10,'LineWidth',1);
plot(Num_targets,Prob_location_CC(1,:),'-.o','Color',parulacustom(1,:),'MarkerSize',10,'LineWidth',1);

plot(Num_targets,Prob_location_LS(2,:),'-square','Color',parulacustom(2,:),'MarkerSize',10,'LineWidth',1);
plot(Num_targets,Prob_location_CC(2,:),'-.square','Color',parulacustom(2,:),'MarkerSize',10,'LineWidth',1);

plot(Num_targets,Prob_location_LS(3,:),'->','Color',parulacustom(3,:),'MarkerSize',10,'LineWidth',1);
plot(Num_targets,Prob_location_CC(3,:),'-.>','Color',parulacustom(3,:),'MarkerSize',10,'LineWidth',1);

plot(Num_targets,Prob_location_LS(4,:),'-x','Color',parulacustom(4,:),'MarkerSize',10,'LineWidth',1);
plot(Num_targets,Prob_location_CC(4,:),'-.x','Color',parulacustom(4,:),'MarkerSize',10,'LineWidth',1);

plot(Num_targets,Prob_location_LS(5,:),'-<','Color',parulacustom(5,:),'MarkerSize',10,'LineWidth',1);
plot(Num_targets,Prob_location_CC(5,:),'-.<','Color',parulacustom(5,:),'MarkerSize',10,'LineWidth',1);

plot(Num_targets,Prob_location_LS(6,:),'-*','Color',parulacustom(6,:),'MarkerSize',10,'LineWidth',1);
plot(Num_targets,Prob_location_CC(6,:),'-.*','Color',parulacustom(6,:),'MarkerSize',10,'LineWidth',1);


legend([lo lx L(1) L(2) L(3) L(4) L(5) L(6)],'Pair-LS','Pair-CC-ESPRIT','$\sigma_e=0$','$\sigma_e=0.1$','$\sigma_e=0.2$','$\sigma_e=0.3$','$\sigma_e=0.4$','$\sigma_e=0.5$','Interpreter','latex','Location','southeast');
xlabel('Number of targets','Interpreter','latex')
ylabel('Probability of successful pairing','Interpreter','latex')

%xlim([2 6])
% xticklabels([2 3 4 5])
% xticks([2 3 4 5])

grid on
set(gca,'GridLineStyle','--')
set(gca,'FontSize',18) 
set(gca,'TickLabelInterpreter','latex')

colormap(parula(length(sigma)))
cb = colorbar;
caxis([0 0.5])
ylabel(cb,'$\sigma_e$','Interpreter','latex','FontSize',18)
%cb.TickLabels={"l",'b','c','d','e'};
%cb.TickLabelInterpreter="latex";

orient tall
print -deps -r300 ProbvsNum_b.eps
