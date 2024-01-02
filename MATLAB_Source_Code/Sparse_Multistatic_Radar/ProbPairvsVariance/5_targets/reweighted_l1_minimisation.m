function [x,solver_status] = reweighted_l1_minimisation(z,P,eb,G,max_iter,epsilon) 
w=ones(G,1);
w_old=w;
for i = 1:max_iter
    W=diag(w);
cvx_begin 
variable xhat(G,1) 
minimize(norm(W*xhat,1));
subject to
norm( z-real(P)*xhat ) <= eb
cvx_end

 w = 1./(epsilon + abs(xhat));
  if norm(w-w_old, 2) < 1e-6
    break
  end
  w_old = w;
end
x = abs(xhat);
%x(abs(x)<1e-8) = 0;
solver_status=cvx_status;
