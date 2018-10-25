clc
clear
%%%%%%%%%%%%%%%%%%产生待估计序列 p(k+1)=a*p(k)+w(t)
a=0.8;   
q=2;          %%w(t)的方差
p0=3;p=[];
for i=1:4000
    p=[p p0];
    p0=a*p0+sqrt(q)*randn(1);
end
 
sum_ar=0;
sum_qr=0;
sum_aY=0;
sum_qY=0;

J=500;
for j=1:J  %%%%%%%%%%%%利用500次估计来计算平均结果
%%%%%%%%%%%%%%%%%%%%%%%%利用最小二乘递推方法得到的估计参数
[ar,qr]=LSRecursive(p);
arlast=ar(max(length(ar)));
qrlast=qr(max(length(qr)));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%利用Yule-Walker方法得到的估计参数
[aY,qY]=YuleWalker(p);
aYlast=aY(max(length(aY)));
qYlast=qY(max(length(qY)));
 
sum_ar=sum_ar+arlast;
sum_qr=sum_qr+qrlast;
sum_aY=sum_aY+aYlast;
sum_qY=sum_qY+qYlast;
end
 
arl=sum_ar/J%%%%求平均值
qrl=sum_qr/J
aYl=sum_aY/J
qYl=sum_qY/J
 
subplot(2,1,1),plot(ar,'k','LineWidth',2)
xlabel('k'),ylabel('方程系数a的估计')
subplot(2,1,2),plot(qr,'k','LineWidth',2)
xlabel('k'),ylabel('方差估计')
figure
subplot(2,1,1),plot(aY,'k','LineWidth',2)
xlabel('k'),ylabel('方程系数a的估计')
subplot(2,1,2),plot(qY,'k','LineWidth',2)
xlabel('k'),ylabel('方差估计')