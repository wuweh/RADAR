clc;
clear;
t=0:0.01:9.99;
m=2*t+20;
v=2;
s3=m+sqrt(v)*randn(1,1000);
plot(t,s3)
xlabel('t')
ylabel('measurement data')
Hk=[t;ones(1,1000)]';
W=0.5*eye(1000);
estim=inv(Hk'*W*Hk)*Hk'*W*s3';
estim