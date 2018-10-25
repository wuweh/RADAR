%转移噪声方差矩阵对滤波效果带来的影响！！！！！

clc;
clear all;
close all;

T = 0.2;
SMP = 50;
F=[ 1 0 T 0; 
    0 1 0 T; 
    0 0 1 0; 
    0 0 0 1]; %状态转移矩阵

dr1 = 0.4;
Q=[ dr1^2 0 0 0; 
    0 dr1^2 0 0; 
    0 0 dr1^2 0; 
    0 0 0 dr1^2]; %转移噪声方差矩阵

H=[ 1 0 0 0; 
    0 1 0 0]; %量测矩阵

dr2 = 2;
R=[ dr2^2 0 ; 
    0 dr2^2]; %量测噪声方差矩阵

P=[ 1 0 0 0; 
    0 1 0 0; 
    0 0 1 0; 
    0 0 0 1]; %协方差矩阵

X(:,1)=[200;200;10;10]; %状态
Z(:,1)= H*X(:,1)+dr2*randn(2,1);
for i=2:SMP  
    X(:,i)= F*X(:,i-1);                          %模拟，产生目标运动真实轨迹
    Z(:,i)= H*X(:,i)+dr2*randn(2,1);     %测量值加上噪声
end

x = X(:,1);
for i=1:SMP
    [x,P] = kf(F,x,P,H,Z(:,i),Q,R);
    x_filter(:,i) = x;
end

 delta_r= sqrt( (x_filter(1,:)-X(1,:)).^2 + (x_filter(2,:)-X(2,:)).^2 );
 
figure;
plot(X(1,:),X(2,:),'-b.');hold on;
plot(x_filter(1,:),x_filter(2,:),'-r.');hold on;
plot(Z(1,:),Z(2,:),'*');hold on;grid on
xlabel('x(m)');ylabel('y(m)');
legend('实际值','滤波值','量测值');

figure;
plot(delta_r,'-b.');hold on;grid on
xlabel('x(m)');
ylabel('y(m)')


