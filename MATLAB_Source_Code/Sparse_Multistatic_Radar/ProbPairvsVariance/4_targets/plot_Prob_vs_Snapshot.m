%% Loading the data
load("Snapshot_test_300.mat")
Prob_location_LS=prob_correct_pairing_LS_1;
Prob_location_CC=prob_correct_pairing_CC_1;


%% Plotting
figure()
subplot(2,1,1)
plot(Q_vector,Prob_location_LS,'--o',Q_vector,Prob_location_CC,'--x','MarkerSize',10)
set ( gca, 'YScale', 'log')

legend('Pair-LS','Pair-CC-ESPRIT','Interpreter','latex') 
legend('Location','southeast')
xlabel('Number of pulses/snapshots','Interpreter','latex')
ylabel('Probability of successful pairing','Interpreter','latex')
xticks={'2/4'; '4/8'; '6/12'; '8/16'};
set(gca,'xtick',[2 4 6 8],'xticklabel',xticks)
%grid on
%set(gca,'GridLineStyle','--')
set(gca,'FontSize',18) 
orient tall
print -deps ProbvsSnapshot.eps
 



