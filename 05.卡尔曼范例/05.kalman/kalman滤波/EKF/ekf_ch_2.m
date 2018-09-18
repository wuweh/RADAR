%扩展Kalman滤波在目标跟踪中的应用
%function EKF_For_TargetTracking
clc;clear;
T=1;%雷达扫描周期

N=60/T;%总的采样次数

X=zeros(4,N);%目标真实位置、速度

X(:,1)=[-100,2,200,20];%目标初始位置、速度
 	
Z=zeros(1,N);%传感器对位置的观测

delta_w=1e-3;%如果增大这个参数，目标的真实轨迹就是曲线了

Q=delta_w*diag([0.5,1]);%过程噪声方差

G=[T^2/2,0;T,0;0,T^2/2;0,T];%过程噪声驱动矩阵
	
R=5;%观测噪声方差
 
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];%状态转移矩阵

x0=200;%观测站的位置

y0=300;
 	
Xstation=[x0,y0];

for t=2:N
    X(:,t)=F*X(:,t-1)+G*sqrtm(Q)*randn(2,1);%目标真实轨迹	
end

for t=1:N	
    Z(t)=Dist(X(:,t),Xstation)+sqrtm(R)*randn;%观测值
end	

%EKF
Xekf=zeros(4,N);	
Xekf(:,1)=X(:,1);%Kalman滤波的状态初始化	
P0=eye(4);%误差协方差矩阵的初始化

for i=2:N
    Xn=F*Xekf(:,i-1);%一步预测
    P1=F*P0*F'+G*Q*G';%预测误差协方差
    dd=Dist(Xn,Xstation);%观测预测	
    %求解雅克比矩阵H	
    H=[(Xn(1,1)-x0)/dd,0,(Xn(3,1)-y0)/dd,0];%泰勒展开的一阶近似
    K=P1*H'*inv(H*P1*H'+R);%卡尔曼最优增益
    Xekf(:,i)=Xn+K*(Z(:,i)-dd);%状态更新
    P0=(eye(4)-K*H)*P1;%滤波误差协方差更新
end
	
%误差分析	
for i=1:N
    Err_KalmanFilter(i)=Dist(X(:,i),Xekf(:,i));%
end 	
%画图	
figure	
hold on;box on;	
plot(X(1,:),X(3,:),'-k');%真实轨迹	
plot(Xekf(1,:),Xekf(3,:),'-r');%扩展Kalman滤波轨迹	
legend(' real trajectory','EKF trajectory');	
xlabel('X-axis  X/m');	
ylabel('Y-axis Y/m');
	
figure 	
hold on ;box on;	
plot(Err_KalmanFilter,'-ks','MarkerFace','r')
xlabel('Time /s');
ylabel('Position estimation bias   /m');

%子函数 求欧氏距离	
function d=Dist(X1,X2);	
    if length(X2)<=2 	
        d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);	
    else 	
        d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(3))^2);
    end
end
	
% %子函数 求欧氏距离
% function d=Dist(X1,X2)
%     if length(X2)<=2
%         d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
%     else	
%         d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(3))^2); 	
%     end
% end