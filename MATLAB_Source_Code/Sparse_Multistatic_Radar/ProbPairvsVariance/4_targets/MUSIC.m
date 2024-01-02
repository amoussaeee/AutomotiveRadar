function [P] = MUSIC(Y,A,Num_Targets,G,C)
%This function returns the MUSIC spectrum using the steering matrix
%provided.
Y=Y(:,1:C);
R=Y*Y';                          %Calculate covariance matrix                         
[Q,D]=eig(R);                    %Eignvalue Decomposition
[~,V]=sort(diag(D),1,'descend'); %Find  largest eigenvalues
Q=Q (:,V);                       %Sort the eigenvectors to put signal eigenvectors first
Qn=Q(:,Num_Targets+1:end);       %Get the noise eigenvectors
noise=(Qn)*Qn';                  %Noise Subspace

P=zeros(G,1); %Initiate the MUSIC spectrum matrix
for g=1:G  %Start loop for parameter search
    matchedsteering=A(:,g);
    P(g,1)=1/(matchedsteering'*noise*matchedsteering); %MUSIC pseudo-spectrum
end                 %End loop for parameter search
end