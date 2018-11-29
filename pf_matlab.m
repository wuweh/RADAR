clc;
clear all;
close all;
x = 0; %初始值
R = 1;
Q = 1;
tf = 100; %跟踪时长
N = 100; %粒子个数
P = 2;
xhatPart = x;
for i = 1 : N 
 xpart(i) = x + sqrt(P) * randn;%初始状态服从0均值，方差为sqrt(P)的高斯分布
end
xArr=[x];
yArr = [x^2 / 20 + sqrt(R) * randn];
xhatArr = [x];
PArr = [P];
xhatPartArr = [xhatPart];
for k = 1 : tf 

x = 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
%k时刻真实值
 y = x^2 / 20 + sqrt(R) * randn; %k时刻观测值
for i = 1 : N
xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) ...
+ 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;%采样获得N个粒子
ypart = xpartminus(i)^2 / 20;%每个粒子对应观测值
vhat = y - ypart;%与真实观测之间的似然
q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R);
%每个粒子的似然即相似度
end
qsum = sum(q);
for i = 1 : N
q(i) = q(i) / qsum;%权值归一化
end
for i = 1 : N %根据权值重新采样
    u = rand;
qtempsum = 0;
 for j = 1 : N

qtempsum = qtempsum + q(j);
 if qtempsum >= u
 xpart(i) = xpartminus(j);
 break;
 end
 end
 end
xhatPart = mean(xpart);
%最后的状态估计值即为N个粒子的平均值，这里经过重新采样后各个粒子的权值相同
xArr = [xArr x]; 
yArr = [yArr y]; 
% xhatArr = [xhatArr xhat];?
PArr = [PArr P];
xhatPartArr = [xhatPartArr xhatPart];
end
t = 0 : tf;
figure;
plot(t, xArr, 'b-.', t, xhatPartArr, 'k-');
legend('Real Value','Estimated Value');
set(gca,'FontSize',10);
xlabel('time step');
ylabel('state');
title('Particle filter')
xhatRMS = sqrt((norm(xArr - xhatArr))^2 / tf);
xhatPartRMS = sqrt((norm(xArr - xhatPartArr))^2 / tf);
figure;
plot(t,abs(xArr-xhatPartArr),'b');
title('The error of PF')