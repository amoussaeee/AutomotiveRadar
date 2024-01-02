function [W_k,F_propagator,F_music] = esprit_propagator_music(Rcc,av_g,cov_length,Num_Targets,fc,W_length)
c=3e8;

%% Propagator
G=Rcc(:,1:Num_Targets);
H=Rcc(:,Num_Targets+1:end);
P_csm=(G'*G)\(G'*H);

Q=[P_csm;-eye(cov_length-Num_Targets)];
Qtotal=Q*Q';
F_propagator=zeros(W_length,1);
for g=1:W_length
    matched_steering=av_g(:,g);
    F_propagator(g)=1/(matched_steering'*Qtotal*matched_steering);
end

%% MUSIC
[Q,D]=eig(Rcc);                    %Eignvalue Decomposition
[~,V]=sort(diag(D),1,'descend'); %Find  largest eigenvalues
Q=Q (:,V);                       %Sort the eigenvectors to put signal eigenvectors first
Qn=Q(:,Num_Targets+1:end);       %Get the noise eigenvectors
noise=(Qn)*Qn';                  %Noise Subspace

F_music=zeros(W_length,1); %Initiate the MUSIC spectrum matrix
for g=1:W_length  %Start loop for parameter search
    matched_steering=av_g(:,g);
    F_music(g,1)=1/(matched_steering'*noise*matched_steering); %MUSIC pseudo-spectrum
end                 %End loop for parameter search



%% ESPRIT
Qs = Q(:,1:Num_Targets); %Signal Subspace

% W1 and W2 matrices:
W1 = Qs(1:cov_length-1, :);
W2 = Qs(2:end, :);

% solve optomization by LS:
psi = inv(W1'*W1) * W1' * W2;

% perform eigendecomposition of psi:
[~, D_psi] = eig(psi);
angles = angle(diag(D_psi));
W_k = angles / (-2*pi*(fc/c));
end

