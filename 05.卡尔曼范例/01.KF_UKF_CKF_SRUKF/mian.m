%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Change Records
%   20180629：1）First upload
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
n=6; 
t=0.05;

global Q R R1 fai gama kesi w m;

m=2*n;
w=1/m;
kesi1 = eye(n);
kesi2 = -eye(n);
kesi = [kesi1,kesi2]*(sqrt(m/2)); 


% Q=[1 0 0 0 0 0;
%     0 1 0 0 0 0;
%     0 0 0.001 0 0 0;
%     0 0 0 1 0 0;
%     0 0 0 0 1 0;
%     0 0 0 0 0 0.001];
% 
% 
R = [2^2    0;
     0      0.002^2];


Q=0.2*eye(6);
    
R1 = [ 2^2    0;
       0    2^2];

f=@(x)[ x(1)+t*x(3)+0.5*t^2*x(5);
        x(2)+t*x(4)+0.5*t^2*x(6); 
        x(3)+t*x(5);
        x(4)+t*x(6);
        x(5);
        x(6)];
    
h=@(x)[ sqrt(x(1)^2+x(2)^2);
        atan(x(2)/x(1))];
    
F = [1 0 t 0 0.5*t*t 0;
     0 1 0 t 0 0.5*t*t;
     0 0 1 0 t 0;
     0 0 0 1 0 t;
     0 0 0 0 1 0;
     0 0 0 0 0 1];
 
 H = [1 0 0 0 0 0;
      0 1 0 0 0 0];

fai = f;
gama = h;

%X,Y,VX,VY,Ax,Ay
s=[70;70;10;20;10;0];

x0=s+sqrtm(Q)*randn(n,1);


P0 = eye(6);
% P0 =[   1 0 0 0 0 0;
%         0 1 0 0 0 0;
%         0 0 0.01 0 0 0;
%         0 0 0 0.01 0 0;
%         0 0 0 0 0.0001 0;
%         0 0 0 0 0 0.0001];

N=50;

Xkf = zeros(n,N);
Xukf = zeros(n,N);
Xckf = zeros(n,N);
X = zeros(n,N);
X1 = zeros(n,N);
X2 = zeros(n,N);
Z = zeros(2,N);
Z1 = zeros(2,N);

for i=1:N
    X(:,i)= f(s);%+sqrtm(Q)*randn(6,1); 
    s = X(:,i);
end

%实际位置显示
% figure
% t=1:N;
% hold on;box on;
% plot(X(1,t),X(2,t),X1(1,t),X1(2,t),X2(1,t),X2(2,t))  

ux=x0;
ux1 = x0;
ux2 = x0;

P1 = P0;
P2 = P0;
%  P0 = chol(P0);
for k=1:N
    Z(:,k)= h(X(:,k)) + sqrtm(R)*randn(2,1);       
    Z1(:,k)= H*X(:,k) + sqrtm(R1)*randn(2,1);       
    
%   [Xukf(:,k),P0] = srukf(f,ux,P0,h,Z(:,k),Q,R);      
    [Xukf(:,k),P0] = ukf(f,ux,P0,h,Z(:,k),Q,R);       
    [Xckf(:,k),P1] = ckf(ux1,P1,Z(:,k));           
%     [Xkf(:,k),P2] = kf(F,ux2,P2,H,Z1(:,k),Q,R1);    
  
    ux  = Xukf(:,k);
    ux1 = Xckf(:,k);
    ux2 = Xkf(:,k);
end

% for k=1:N
%     RMS_R = sqrt(X(1,k)^2+X(2,k)^2);
%     RMS(k)  =   RMS_R - sqrt(Xukf(1,k)^2+Xukf(2,k)^2);
%     RMS1(k) =   RMS_R - sqrt(Xckf(1,k)^2+Xckf(2,k)^2);
%     RMS2(k) =   RMS_R - sqrt(Xkf(1,k)^2+Xckf(2,k)^2);
%     RMS3(k) =   RMS_R - Z(1,k);
%     RMS4(k) =   RMS_R - sqrt(Z1(1,k)^2+Z1(2,k)^2);    
%     
%     RMS_A = atan(X(1,k)/X(2,k));
%     RMS_1(k)    =   RMS_A - atan(Xukf(1,k)/Xukf(2,k));
%     RMS1_1(k)   =   RMS_A - atan(Xckf(1,k)/Xckf(2,k));
%     RMS2_1(k)   =   RMS_A - atan(Xkf(1,k)/Xkf(2,k));
%     RMS3_1(k)   =   RMS_A - Z(2,k);
%     RMS4_1(k)   =   RMS_A - atan(Z1(1,k)/Z1(2,k));
% end

% figure
% t=1:N;
% hold on;
% box on;
% plot(X(1,t),X(2,t),'-k.')  
% plot(Z(1,t).*cos(Z(2,t)),Z(1,t).*sin(Z(2,t)),'-b.')
% % plot(Z1(1,t),Z1(2,t),'-r.')
% % axis([0 15 0 200]);
% xlabel('x/m');
% ylabel('y/m');
% title('Real And Measure Position');
% legend('Real Position','Measure 1','Measure 2');


figure
t=1:N;
hold on;box on;
plot(X(1,t),X(2,t),'-b.')        
plot(Xukf(1,t),Xukf(2,t),'-r.');
plot(Z(1,t).*cos(Z(2,t)),Z(1,t).*sin(Z(2,t)),'*')
plot(Xckf(1,t),Xckf(2,t),'-b.');
% plot(Xkf(1,t),Xkf(2,t),'-g.');
xlabel('x/m');
ylabel('y/m');
% legend('Target','UKF','CKF','KF');
title('UKF');




