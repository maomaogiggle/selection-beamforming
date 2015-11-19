function [v,u,g] = opdual(W,N,K)
cvx_begin
    variable u(N); variable v; variable g;  
    maximize g;
    subject to
       [W + diag(u)+ v*eye(N),-0.5*u;-0.5*u',-K*v-g] == semidefinite(N+1);
cvx_end
end