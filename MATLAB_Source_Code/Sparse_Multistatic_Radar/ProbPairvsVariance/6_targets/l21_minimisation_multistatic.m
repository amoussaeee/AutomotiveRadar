function [estimated,solver_status] = l21_minimisation_multistatic(Z,A_1,A_2,eb,G,P)
cvx_begin 
variable X_1(G,P)
variable X_2(G,P)
expression x(G,1)

for index = 1:G
            u = [X_1(index,:),X_2(index,:)];
            x(index,1) = norm(u,2);
end

minimize(sum(x));
subject to
norm( Z-[real(A_1)*X_1, real(A_2)*X_2] ,'fro' ) <= eb
cvx_end
estimated=x;
solver_status=cvx_status;
end
