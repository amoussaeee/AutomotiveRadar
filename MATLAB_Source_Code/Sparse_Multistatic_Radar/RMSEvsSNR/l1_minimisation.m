function [estimated,solver_status] = l1_minimisation(z,P,eb,G) 
cvx_begin 
variable x(G,1) 
minimize(norm(x,1));
subject to
norm( z-real(P)*x ) <= eb
cvx_end
estimated=abs(x);
solver_status=cvx_status;
end