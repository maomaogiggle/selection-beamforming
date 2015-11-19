% ==================================================================== %
% This code is used to get the trade-off figure 6
% Author: Xiangrong
% Date: 25/11/2012
% ==================================================================== %
clear;
clc;
close all;

N = 4;
K = 2:N*N;
n = -(N-1)/2:1:(N-1)/2;
m = [-(N-1)/2:1:(N-1)/2]';
mn = [];
for i = 1:length(n)
    mn = [mn;[n(i)*ones(N,1),m]]; %mn is the position vector of N*N/2 by 2 dimension
end
lambda = 1;
d = lambda/2;
k0 = 2*pi/lambda;
phis = 0.2*pi;
thetas = 0.2*pi;
us = [cos(thetas)*cos(phis);cos(thetas)*sin(phis)];
%computational cost
cost = (K.^3)/((N*N)^3);

thetaj = 0.22*pi;
phij = 0.22*pi;
uj = [cos(thetaj)*cos(phij);cos(thetaj)*sin(phij)];
y1 = zeros(length(K),1);
SINR1 = zeros(length(K),1);
for k = 1:length(K)
 w = exp(1i*k0*d*mn*(us-uj));
 W = w*w';
 W = (1/(K(k)^2))*W;
 [y1(k),x0,X0] = sdprelaxation(W,N*N,K(k));
 SINR1(k) = K(k)*(1-y1(k));
end

SINR_n1 = SINR1./max(SINR1);

%scenario 2
thetaj = 0.4*pi;
phij = 0.4*pi;
uj = [cos(thetaj)*cos(phij);cos(thetaj)*sin(phij)];
y2 = zeros(length(K),1);
SINR2 = zeros(length(K),1);
for k = 1:length(K)
 w = exp(1i*k0*d*mn*(us-uj));
 W = w*w';
 W = (1/(K(k)^2))*W;
 [y2(k),x0,X0] = sdprelaxation(W,N*N,K(k));
 SINR2(k) = K(k)*(1-y2(k));
end

SINR_n2 = SINR2./max(SINR2);

%plotting
plotyy(cost,10*log10(SINR_n1),cost,10*log10(SINR_n2));

