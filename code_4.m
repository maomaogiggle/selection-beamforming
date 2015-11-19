%testing CM method
clear;
clc;
close all;

N = 4;
K = 8;
n = -(N-1)/2:1:(N-1)/2;
m = [-(N-1)/2:1:(N-1)/2]';
mn = [];
for i = 1:length(n)
    mn = [mn;[n(i)*ones(N,1),m]]; %mn is the position vector of N*N/2 by 2 dimension
end
index = 1:N*N; % index is the index for mn
lambda = 1;
d = lambda/2;
k0 = 2*pi/lambda;
phis = 0.3*pi;
thetas = 0.2*pi;
us = [cos(thetas)*cos(phis);cos(thetas)*sin(phis)];
thetaj = 0.3*pi;
phij = [0:0.01:0.5]*pi;
uj = [cos(thetaj)*cos(phij);cos(thetaj)*sin(phij)];

V = nchoosek(index,K);
num = nchoosek(N*N,K);
e = ones(1,K);
minscc = zeros(length(phij),1);
y = zeros(length(phij),1);
alpha = zeros(length(phij),1);

for k = 1:length(phij)
    scc = zeros(num,1);
    for i = 1:num
        Po = [];
        for q = 1:K
            Po = [Po;mn(V(i,q),:)];
        end
        vsj = exp(1i*k0*d*Po*(us-uj(:,k)));
        scc(i) = abs(e*vsj/K)^2;
    end
 minscc(k) = min(scc);
 w = exp(1i*k0*d*mn*(us-uj(:,k)));
 W = w*w';
 W = (1/(K^2))*W;
 [kesi_v_opt,kesi_mu_opt,y(k)] = opdual(real(W),N*N,K);
 [alpha(k),x0] = CM(real(W),K,N*N);
end

figure;
plot(phij,minscc,'ro-','LineWidth',2);
hold on;
plot(phij,y,'b+-','LineWidth',2);
hold on;
plot(phij,alpha,'kh:','LineWidth',2);
xlabel('arrival direction of interference: rad');
ylabel('squared scc value');
legend('true minimum scc','lower bound from Lagrange Dual','CM');


figure;
plot(phij,alpha-minscc);
hold on;
plot(phij,1/(K*K)*ones(1,length(phij)),'r');
legend('distance between CM and the true','upper bound');


 