clc
clear
%%%测量参数：周期和测量方差的设置
T=1;
R=40^2;;
%%%圆形机动的测量值
Tt=31;Tz=57;[t,y]=funtrackinglinecircle(T,Tt,Tz,R);
%%%%%%%%%%选择模型
qq=1;%%%此值的选择对性能影响很大，可以去不同的值，如取0.1，1，10等，直接修改程序中的值即可。

 [A,Q]=CVmodel(T,qq);C=[1 0]; 
% [A,Q]=CAmodel(T,qq);C=[1 0 0]; %%这两个模型在执行的时候选其中的一个。
%使用Kalman进行滤波,对于圆形机动，需要横纵轴分别估计
%%%估计横轴
xe=zeros(length(Q),1);p=1000*eye(size(A));xx1=[];
for i=1:length(t)
[xe,p]=kalmanfun(A,C,Q,R,xe,y(1,i),p);
xx1=[xx1 xe];
end
%%%%估计纵轴
xe=zeros(length(Q),1);p=1000*eye(size(A));xx2=[];
for i=1:length(t)
[xe,p]=kalmanfun(A,C,Q,R,xe,y(2,i),p)
xx2=[xx2 xe];
end
plot(y(1,:),y(2,:),'*');hold on,plot(C*xx1,C*xx2);hold off
figure