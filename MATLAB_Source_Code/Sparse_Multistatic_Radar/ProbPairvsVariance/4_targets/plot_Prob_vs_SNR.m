%% Loading the data
load("SNR_test_500.mat")
Prob_location_LS=prob_correct_pairing_LS_1;
Prob_location_CC=prob_correct_pairing_CC_1;


%% Plotting
figure()
subplot(2,1,1)
plot(SNR_vector,Prob_location_LS,'--o',SNR_vector,Prob_location_CC,'--x','MarkerSize',10)


set ( gca, 'YScale', 'log')

legend('Pair-LS','Pair-CC-ESPRIT','Interpreter','latex') 
legend('Location','southeast')
xlabel('SNR [dB]','Interpreter','latex')
ylabel('Probability of successful sairing','Interpreter','latex')
xlim([120 160])
yticklabels([0.01 0.1 1])
%grid on
%set(gca,'GridLineStyle','--')
set(gca,'FontSize',18) 
orient tall
print -deps ProbvsSNR.eps
 



