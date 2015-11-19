%=========================================================================%
% This code is CM heuristic method (Correlation Measurement)
% Author: Xiangrong Wang
% Date:   20/10/2012
%=========================================================================%
function [alpha,x] = CM(W,K,N)
temp_W = W;
x = ones(N,1);
for k = 1:N-K
    e = ones(N,1);
    cm = zeros(N,1);
    for i = 1:N
        cm(i) = e'*temp_W(:,i);
    end
    [Y,I] = max(cm);
    temp_W(:,I) = 0;
    temp_W(I,:) = 0;
    x(I) = 0;
end
W_c = temp_W;
alpha = ones(1,N)*W_c*ones(N,1);
end
















