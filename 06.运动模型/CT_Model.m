clc;
close all;
clear all;

T = 0.1;
w = 20/180*pi;  %转向角速度,单位为弧度
SPM = 100;
dx = 1;
dy = 1;
dvx = 0.5;
dvy = 0.5;

%匀转角运动模型
F = [    1,       sin(w*T)/w,           0,                  -1*(1-cos(w*T))/w,      0;
        0,      cos(w*T),               0,                  -1*sin(w*T),                0;
        0,      (1-cos(w*T))/w,      1,                   sin(w*T)/w,               0;
        0,      sin(w*T),                0,                   cos(w*T),                  0;
        0,      0,                          0,                   0,                             1];

%匀速直线运动模型
F1 = [ 1,  T,  0,  0;
      0,  1,  0,  0;
      0   0,  1,  T;
      0 , 0 , 0 , 1;];

H = [1,0,0,0,0;
     0,0,1,0,0 ];

H1 = [1,0,0,0;
      0,0,1,0 ];

P2 =[ 1 0 0 0 0 ;
      0 1 0 0 0 ;
      0 0 1 0 0 ;
      0 0 0 1 0 ;
      0 0 0 0 1;  ];

P21 =[  1 0 0 0 ;
        0 1 0 0 ;
        0 0 1 0  ;
        0 0 0 1  ; ];

Q1 = [ dx^2,0,0,0;
        0,dvx^2,0,0;
        0,0,dy^2,0;
        0,0,0,dvy^2; ];

Q = [dx^2,0,0,0,0;
     0,dvx^2,0,0,0;
     0,0,dy^2,0,0;
     0,0,0,dvy^2,0;
     0,0,0,0,0.01; ];

R = [dx^2   0;
        0    dy^2];

R1 = R;
  

%初始化
x(:,1) = [20,10,20,10,w]';
Z1(:,1) = H*x(:,1)+sqrtm(R)*rand(2,1);

%生成运动轨迹
for i = 2:SPM
    x(:,i) = F*x(:,i-1);%+sqrtm(Q)*rand(5,1);
    Z1(:,i) = H*x(:,i)+sqrtm(R)*rand(2,1);
end

ux2 = x(:,1);
ux3 = x(1:4,1);
for k = 1:SPM
    [Xkf(:,k),P2] = kf(F,ux2,P2,H,Z1(:,k),Q,R);            %调用kf算法 
    [Xkf1(:,k),P21] = kf(F1,ux3,P21,H1,Z1(:,k),Q1,R1);     %调用kf算法 
    ux2 = Xkf(:,k);
    ux3 = Xkf1(:,k);
end
 
figure;
subplot(221)
hold on;
plot(x(1,:),x(3,:),'-b.');
hold on;
plot(Z1(1,:),Z1(2,:),'g*');
plot(Xkf(1,:),Xkf(3,:),'-r.');
plot(Xkf1(1,:),Xkf1(3,:),'-y.');
% set(gca,'color',[0.3, 0.3, 0.3]);
legend('真实值','测量值','CT预测值','CV预测值')
  
D=sqrt((x(1,:)-Xkf(1,:)).^2+(x(3,:)-Xkf(3,:)).^2);%均方根误差
D1=sqrt((x(1,:)-Xkf1(1,:)).^2+(x(3,:)-Xkf1(3,:)).^2);%均方根误差
subplot(222)
plot(D(1:SPM),'-bo','MarkerSize',4);
hold on
plot(D1(1:SPM),'-ro','MarkerSize',4);
legend('CT模型误差','CV模型误差')

x_p = 1:1:SPM;
delta_x1=x(1,:)-Xkf(1,:);
delta_y1=x(3,:)-Xkf(3,:);
delta_x2=x(1,:)-Xkf1(1,:);
delta_y2=x(3,:)-Xkf1(3,:);

subplot(223)
plot(x_p,delta_x1,'-ro','MarkerSize',4);
hold on
plot(x_p,delta_x2,'-bo','MarkerSize',4);

subplot(224)
plot(x_p,delta_y1,'-ro','MarkerSize',4);
hold on
plot(x_p,delta_y2,'-bo','MarkerSize',4);
 
 
