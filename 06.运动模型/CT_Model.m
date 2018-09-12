clc;
close all;
clear all;

T = 0.5;
w = -20/180*pi;  %转向角速度,单位为弧度
SPM = 100;

F = [1,sin(w*T)/w, 0, -1*(1-cos(w*T))/w,0;
     0,cos(w*T),0,-1*sin(w*T), 0;
     0,(1-cos(w*T))/w,1,sin(w*T)/w,0;
     0,sin(w*T),0,cos(w*T),0;
     0,0,0,0,1];
 
 F1 = [1,T,0,0;
        0,1,0,0;
        0,0,1,T;
        0 ,0 ,0 ,1;];
 
 H = [1,0,0,0,0;
      0,0,1,0,0];
   H1 = [1,0,0,0;
      0,0,1,0];
 
 P2 =[  1 0 0 0 0 ;
        0 1 0 0 0 ;
        0 0 1 0 0 ;
        0 0 0 1 0 ;
        0 0 0 0 1;
        ];
 
 Q = [4,0,0,0,0;
      0,1,0,0,0;
      0,0,4,0,0;
      0,0,0,1,0;
      0,0,0,0,0.01;
     ];
 
 R = [3^2 0;
      0   3^2];
  
   P21 =[1 0 0 0 ;
        0 1 0 0 ;
        0 0 1 0  ;
        0 0 0 1  ;
        ];
 
 Q1 = [4,0,0,0;
      0,1,0,0;
      0,0,4,0;
      0,0,0,1;
     ];
 
 R1 = [3^2 0;
      0   3^2];
  
 %初始化
 x(:,1) = [20,10,20,10,w]';
 x2(:,1) = [20,10,20,10];
 Z1(:,1) = H*x(:,1)+sqrtm(R)*rand(2,1);
 ux2 = x(:,1);
 ux3 = x2(:,1);
 
 %
 for i = 2:SPM
    x(:,i) = F*x(:,i-1)+sqrtm(Q)*rand(5,1);
    Z1(:,i) = H*x(:,i)+sqrtm(R)*rand(2,1);
    x2(1:4,i) = x(1:4,i);
 end
 
 figure;
 hold on;
 plot(x(1,:),x(3,:),'-b.');
 hold on;
 plot(Z1(1,:),Z1(2,:),'-r.');
 
sumx1 = 0;
sumx2 = 0;
 for k = 1:SPM
     [Xkf(:,k),P2] = kf(F,ux2,P2,H,Z1(:,k),Q,R);     %调用kf算法 
     sumx1 = sumx1 + Xkf(1,k);
     ux2 = Xkf(:,k);
%      ux2(1) = sumx1/k;
     [Xkf1(:,k),P21] = kf(F1,ux3,P21,H1,Z1(:,k),Q1,R1);     %调用kf算法 
     sumx2 = sumx2 + Xkf1(1,k);
     ux3 = Xkf1(:,k);
%      ux3(1) = sumx2/k;

 end
 
plot(Xkf(1,:),Xkf(3,:),'-g.');
plot(Xkf1(1,:),Xkf1(3,:),'-c.');
legend('Real Position','Measured Position','Kalman Prediction')
  
% D=sqrt((x(1,:)-Xkf(1,:)).^2+(x(3,:)-Xkf(3,:)).^2);%均方根误差
% figure
% plot(D(1:SPM),'-r*');

 
 
 
