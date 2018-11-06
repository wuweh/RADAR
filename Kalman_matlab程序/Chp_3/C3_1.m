clc;
clear;
t=0:0.01:9.99;
m=2.5*t+22;
v=2;
s3=m+sqrt(v)*randn(1,1000);
plot(t,s3)
xlabel('t')
ylabel('measurement data')
Hk=[t;ones(1,1000)]';
estim=inv(Hk'*Hk)*Hk'*s3';
estim 
