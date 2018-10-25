clc
clear
A=1;C=1;Q=9;R=4;I=10;p0=10;
[pp1,pp,KK]=steadycov(A,C,Q,R,I,p0);
subplot(3,1,1);plot(pp1);xlabel('k'),ylabel('向前一步预测估计方差')
subplot(3,1,2);plot(pp);xlabel('k'),ylabel('状态估计方差')
subplot(3,1,3);plot(KK);xlabel('k'),ylabel('滤波器增益')