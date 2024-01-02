function [estimated,solver_status] = l21_minimisation(Z,A,eb,G,P)
cvx_begin 
variable X(G,P)
expression x(G,1)

for index = 1:G
            u = [X(index,:)];
            x(index,1) = norm(u,2);
end

minimize(sum(x));
subject to
norm( Z-real(A)*X ,'fro' ) <= eb
cvx_end
estimated=x;
solver_status=cvx_status;
end