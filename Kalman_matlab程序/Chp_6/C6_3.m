clc
clear
x0=[2000,0];
v=300;%线速度
T=1;
r=2000;%半径
w=v/r;
op=[0 0];%圆心
x=[];xt=[];www=[];ww=0;
for t=1:T:50
ww=ww+w;
    x2=[cos(ww)*r+op(1);sin(ww)*r+op(2)];
    x=[x x2];
    xt=[xt t];
    www=[www ww];
end
plot(x(1,:),x(2,:),'*')
figure
plot(xt,x,'*')