function [alpha,x0,X0] = sdprelaxation(W,N,K)
cvx_begin sdp
    variable x(N); variable X(N,N) symmetric;
    alpha = trace(W*X);
    minimize (alpha);
    subject to
       0 <= trace(X)-K <= 0;
       [X,x;x',1] == semidefinite(N+1);
       for i = 1:N
           e = [zeros(i-1,1);ones(1,1);zeros(N-i,1)];
           E = e*e';
           0 <= trace(X*E)-e'*x <= 0;
       end
cvx_end
x0 = x;     %return optimal selection vector
X0 = X;
end