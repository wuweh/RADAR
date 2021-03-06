%****************************************************************************************************%
%比较CT、CV模型
%20181201：角速度越大，CT模型与CV模型的跟踪性能相差约大
%20181212: 当真实角速度与模型角速度不一致时，CT模型的跟踪效果会下降甚至不如CV模型
%****************************************************************************************************%

clc;
close all;
clear all;

T = 0.2;
SPM = 100;
w = 6/180*pi; 
w1 = 6/180*pi;  %转向角速度,单位为弧度  若w接近于0，轨迹近似于直线运动

%匀转角运动模型
F = [1,sin(w*T)/w,0,-1*(1-cos(w*T))/w,0;
         0,cos(w*T),0,-1*sin(w*T),0;
         0,(1-cos(w*T))/w,1,sin(w*T)/w,0;
         0,sin(w*T),0,cos(w*T),0;
         0,0,0,0,1];
     
F2 = [1,sin(w1*T)/w1,0,-1*(1-cos(w1*T))/w1,0;
         0,cos(w1*T),0,-1*sin(w1*T),0;
         0,(1-cos(w1*T))/w1,1,sin(w1*T)/w1,0;
         0,sin(w1*T),0,cos(w1*T),0;
         0,0,0,0,1];

%匀速直线运动模型
F1 = [ 1,  T,  0,  0;
      0,  1,  0,  0;
      0   0,  1,  T;
      0 , 0 , 0 , 1;];
  
H = [1,0,0,0,0; 0,0,1,0,0 ];
H1 = [1,0,0,0; 0,0,1,0 ];

P2 =[ 1 0 0 0 0 ;
      0 1 0 0 0 ;
      0 0 1 0 0 ;
      0 0 0 1 0 ;
      0 0 0 0 1;  ];

P21 =[  1 0 0 0 ;
        0 1 0 0 ;
        0 0 1 0  ;
        0 0 0 1  ; ];

q = 1;
Q_CV=[ T^3/3 T^2/2 0 0; 
       T^2/2 T 0 0; 
        0 0 T^3/3 T^2/2; 
        0 0 T^2/2 T]; %转移噪声方差矩阵
    
Q_CV = Q_CV*q;

qw = 0.1;
Q_CT=[ q*T^3/3 q*T^2/2 0 0 0; 
            q*T^2/2 q*T 0 0 0; 
            0 0 q*T^3/3 q*T^2/2 0; 
            0 0 q*T^2/2 q*T 0;
            0 0 0 0 qw*T]; %转移噪声方差矩阵
    
Q_CA=[ T^5/20 T^4/8 T^3/6 0 0 0; 
           T^4/8 T^3/3 T^2/2 0 0 0;
           T^3/6 T^2/2 T 0 0 0;
            0 0 0 T^5/20 T^4/8 T^3/6; 
            0 0 0 T^4/8 T^3/3 T^2/2;
            0 0 0 T^3/6 T^2/2 T]; %转移噪声方差矩阵

dx = 1;
dy = 1;
R = [dx^2,0;0,dy^2];
R1 = R;
  
%初始化
x(:,1) = [20,10,20,10,w1]';
%生成运动轨迹
for i = 2:SPM
    x(:,i) = F2*x(:,i-1);%+sqrtm(Q)*rand(5,1);
end

%生成量测目标
for i=1:SPM
    Z1(:,i) = H*x(:,i)+sqrtm(R)*randn(2,1);
end

Xkf(:,1) = x(:,1);
Xkf1(:,1) = x(1:4,1);
for k = 1:SPM-1
    [Xkf(:,k+1),P2] = kf(F,Xkf(:,k),P2,H,Z1(:,k+1),Q_CT,R);            %CT
    [Xkf1(:,k+1),P21] = kf(F1,Xkf1(:,k),P21,H1,Z1(:,k+1),Q_CV,R1);     %CV
end
 
figure;
subplot(221)
hold on;
plot(x(1,:),x(3,:),'-b.');hold on;grid on;
plot(Z1(1,:),Z1(2,:),'.');hold on;grid on;
plot(Xkf(1,:),Xkf(3,:),'-r.');hold on;grid on;
plot(Xkf1(1,:),Xkf1(3,:),'-g.');hold on;grid on;
legend('真实值','测量值','CT预测值','CV预测值');
hold on;grid on;
  
D=sqrt((x(1,:)-Xkf(1,:)).^2+(x(3,:)-Xkf(3,:)).^2);%均方根误差
D1=sqrt((x(1,:)-Xkf1(1,:)).^2+(x(3,:)-Xkf1(3,:)).^2);%均方根误差
subplot(222)
plot(D(1:SPM),'-ro','MarkerSize',4);hold on;grid on;
plot(D1(1:SPM),'-bo','MarkerSize',4);hold on;grid on;
legend('CT模型误差','CV模型误差')

x_p = 1:1:SPM;
delta_x1=x(1,:)-Xkf(1,:);
delta_y1=x(3,:)-Xkf(3,:);
delta_x2=x(1,:)-Xkf1(1,:);
delta_y2=x(3,:)-Xkf1(3,:);

subplot(223)
plot(x_p,abs(delta_x1),'-ro','MarkerSize',4);hold on;grid on;
plot(x_p,abs(delta_x2),'-bo','MarkerSize',4);hold on;grid on;
legend('CT模型误差','CV模型误差')

subplot(224)
plot(x_p,abs(delta_y1),'-ro','MarkerSize',4);hold on;grid on;
plot(x_p,abs(delta_y2),'-bo','MarkerSize',4);hold on;grid on;
legend('CT模型误差','CV模型误差')
 
 
