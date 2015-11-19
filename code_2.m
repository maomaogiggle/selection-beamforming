% ==================================================================== %
% This code is used to produce Fig. 5 to investigate the SCC with different K
%
% Author: Xiangrong
% Date: 23/11/2012
% ==================================================================== %
clear;
clc;
close all;

N = 16;
n = [0:N-1]';
K = 2:16;
lambda = 1;
d = lambda/2;
k0 = 2*pi/lambda;
phis = 0.2*pi;
us = cos(phis);
phij = 0.16*pi; %change it to 0.16pi, 0.23pi, 0.25pi
uj = cos(phij);
y = zeros(length(uj),length(K));

for i = 1:length(K)
    for q = 1:length(uj)
        w = exp(1i*k0*d*n*(us-uj(q)));
        W = w*w';
        W = (1/(K(i)^2))*W;
        [y(q,i),x0,X0] = sdprelaxation(W,N,K(i));
    end
end
 
figure;
plot(K,y,'kd-.','LineWidth',2);
hold on;
xlabel('arrival direction of interference: rad');
ylabel('squared scc value');
