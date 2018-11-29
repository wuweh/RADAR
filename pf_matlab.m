clc;
clear all;
close all;
x = 0; %��ʼֵ
R = 1;
Q = 1;
tf = 100; %����ʱ��
N = 100; %���Ӹ���
P = 2;
xhatPart = x;
for i = 1 : N 
 xpart(i) = x + sqrt(P) * randn;%��ʼ״̬����0��ֵ������Ϊsqrt(P)�ĸ�˹�ֲ�
end
xArr=[x];
yArr = [x^2 / 20 + sqrt(R) * randn];
xhatArr = [x];
PArr = [P];
xhatPartArr = [xhatPart];
for k = 1 : tf 

x = 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
%kʱ����ʵֵ
 y = x^2 / 20 + sqrt(R) * randn; %kʱ�̹۲�ֵ
for i = 1 : N
xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) ...
+ 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;%�������N������
ypart = xpartminus(i)^2 / 20;%ÿ�����Ӷ�Ӧ�۲�ֵ
vhat = y - ypart;%����ʵ�۲�֮�����Ȼ
q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R);
%ÿ�����ӵ���Ȼ�����ƶ�
end
qsum = sum(q);
for i = 1 : N
q(i) = q(i) / qsum;%Ȩֵ��һ��
end
for i = 1 : N %����Ȩֵ���²���
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
%����״̬����ֵ��ΪN�����ӵ�ƽ��ֵ�����ﾭ�����²�����������ӵ�Ȩֵ��ͬ
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