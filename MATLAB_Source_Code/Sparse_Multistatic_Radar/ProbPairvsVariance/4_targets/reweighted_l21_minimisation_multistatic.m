function [estimated,solver_status] = reweighted_l21_minimisation_multistatic(Z,A_1,A_2,eb,G,P,max_iter,epsilon)
w=ones(G,1);
w_old=w;
for i = 1:max_iter
    W=diag(w);

cvx_begin 
variable X_1(G,P)
variable X_2(G,P)
expression x(G,1)

for index = 1:G
            u = [X_1(index,:),X_2(index,:)];
            x(index,1) = norm(u,2);
end

minimize(sum(W*x));
subject to
norm( Z-[real(A_1)*X_1, real(A_2)*X_2] ,'fro' ) <= eb
cvx_end

 w = 1./(epsilon + abs(x));
  if norm(w-w_old, 2) < 1e-6
    break
  end
  w_old = w;
end
x(abs(x)<1e-8) = 0;
estimated=x;

solver_status=cvx_status;
end
