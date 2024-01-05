function [estimated,solver_status] = l21_minimisation(Z,A,eb,G,P)
% This function is for solving an l21 minimisation problem using CVX.
% Author: Ali Moussa - Last modified: 05/01/2024

cvx_begin 
variable X(G,P)      
expression x(G,1)

for index = 1:G
            u = [X(index,:)];
            x(index,1) = norm(u,2);
end

minimize(sum(x));
subject to
norm( Z(:,1:P)-real(A)*X ,'fro' ) <= eb
cvx_end
estimated=x;
solver_status=cvx_status;
end