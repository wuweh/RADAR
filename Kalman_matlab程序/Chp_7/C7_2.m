clc
clear
%%%测量参数：周期和测量方差的设置
%使用kalman滤波器对蛇形机动进行估计
%%%产生蛇形机动的测量值
T=0.3;
t=1:T:30;a=10;omig=pi/3;R=2;
[yreal,ym]=funtrackingsnake(a,omig,t,R);
%%%估计
qq=10;%%%此值的选择对性能影响很大，可以去不同的值，如取0.1，1，10等，直接修改程序中的值即可。

[A,Q]=CAmodel(T,qq);C=[1 0 0];

xe=zeros(length(Q),1);p=1000*eye(size(A));xx1=[];
for i=1:length(t)
[xe,p]=kalmanfun(A,C,Q,R,xe,ym(i),p)
xx1=[xx1 xe];
end
plot(t,ym,t,C*xx1)

%%%%%%%%%%选择模型
