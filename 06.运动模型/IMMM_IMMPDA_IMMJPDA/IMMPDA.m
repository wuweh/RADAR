% 基于IMM算法的目标跟踪
%   20181223: 1）更新IMMPDA算法
%   20181225: 1）更新IMMPDA与IMM算法的误差比较
%                    2）注：把三个模型的预测值进行平均处理，会失去IMM的意义，目前表现的和CT模型结果类似
clc;clear all;close all;

tic
n=50;
global Pd Pg gamma G Q R noise_total;
Pd=1;       
Pg=0.97;   
gamma = 1; 

g_sigma = 1;     
  
simTime=100;      %仿真迭代次数
T=0.5;                     %采样时间
w2=5*2*pi/360;     %模型2转弯率3度
w3=-5*2*pi/360;    %模型3转弯率-3度
H=[1,0,0,0;0,0,1,0];                      %模型量测矩阵
G=[T^2/2,0;T,0;0,T^2/2;0,T];              %模型过程噪声加权矩阵
Q=[2^2,0;0,2^2];                                  %模型过程噪声协方差矩阵
r=3^2;                                 %20 2000
R=[r,0;0,r];                            %模型量测噪声协方差矩阵

F1=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];     %模型1状态转移矩阵

F2=[1,sin(w2*T)/w2,0,(cos(w2*T)-1)/w2;
    0,cos(w2*T),0,sin(w2*T);
    0,(1-cos(w2*T))/w2,1,sin(w2*T)/w2;
    0,-sin(w2*T),0,cos(w2*T)];            %模型2状态转移矩阵 左转弯

F3=[1,sin(w3*T)/w3,0,(cos(w3*T)-1)/w3;
    0,cos(w3*T),0,sin(w3*T);
    0,(1-cos(w3*T))/w3,1,sin(w3*T)/w3;
    0,-sin(w3*T),0,cos(w3*T)];            %模型3状态转移矩阵 右转弯

x0=[100,20,100,20]';  % 初始状态

% 产生量测数据
x = zeros(4,simTime);
z = zeros(2,simTime);         

x(:,1)=x0;
z(:,1)=H*x(:,1)+sqrt(R)*randn(2,1);


noise_total = 10; 
z_noise = zeros(noise_total+1,2,simTime);   
ellipse_Volume= pi*gamma*20;
side=sqrt((ellipse_Volume*gamma+1)/gamma)/2;
Noise_x= z(1,1)+side-2*rand(1,noise_total)*side; 
Noise_y= z(2,1)+side-2*rand(1,noise_total)*side;  
z_noise(:,:,1) =[[Noise_x;Noise_y]'; z(:,1)'];

for a=2:simTime
    if (a>=20)&&(a<=50) 
        x(:,a)=F2*x(:,a-1);      
    elseif (a>=50)&&(a<=80) 
        x(:,a)=F3*x(:,a-1);        
    else
        x(:,a)=F1*x(:,a-1);     
    end
    z(:,a)=H*x(:,a)+sqrt(R)*randn(2,1);
    
    Noise_x= z(1,a)+side-2*rand(1,noise_total)*side; 
    Noise_y= z(2,a)+side-2*rand(1,noise_total)*side;  
    z_noise(:,:,a) = [[Noise_x;Noise_y]'; z(:,a)'];
end

%作图部分
figure
plot(x(1,:),x(3,:),'b*-'); hold on;grid on;
for a=1:simTime
    plot(z_noise(:,1,a), z_noise(:,2,a),'k.');
end

x_pro_IMM(:,1)=x0;
x_pro_IMM_PDA(:,1)=x0;
%模型转移概率矩阵
pij=[0.95,0.025,0.025;
       0.025,0.95,0.025;
       0.025,0.025,0.95];    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%基本IMM算法参数
u_IMM=zeros(3,simTime);
%IMM算法模型概率
u_IMM(:,1)=[0.8,0.1,0.1]';  
%IMM算法各模型初始状态
x1_IMM=x0;x2_IMM=x0;x3_IMM=x0; 

 %初始状态协方差矩阵
P0=diag([10,10,10,10]); 
P1_IMM=P0;P2_IMM=P0;P3_IMM=P0;
P_IMM(:,:,1)=P0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMMPDA
u_IMM_PDA=zeros(3,simTime);
u_IMM_PDA(:,1)=[0.8,0.1,0.1]';  
x1_IMM_PDA=x0;x2_IMM_PDA=x0;x3_IMM_PDA=x0; 

P0=diag([10,10,10,10]);  
P1_IMM_PDA=P0;P2_IMM_PDA=P0;P3_IMM_PDA=P0;
P_IMM_PDA(:,:,1)=P0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
x4 = x1_IMM;
x4_rec(:,1) = x4;
P4 = P0;

for t=1:simTime-1
    %计算t时刻处于每个模型的概率
    c_j=pij'*u_IMM(:,t); 
    
    %计算每个模型的条件概率
    ui1=(1/c_j(1))*pij(:,1).*u_IMM(:,t);
    ui2=(1/c_j(2))*pij(:,2).*u_IMM(:,t);
    ui3=(1/c_j(3))*pij(:,3).*u_IMM(:,t);    
    
    % 计算各模型滤波初始化条件
    x01=x1_IMM*ui1(1)+x2_IMM*ui1(2)+x3_IMM*ui1(3);
    x02=x1_IMM*ui2(1)+x2_IMM*ui2(2)+x3_IMM*ui2(3);
    x03=x1_IMM*ui3(1)+x2_IMM*ui3(2)+x3_IMM*ui3(3);   
    
    %各模型滤波初始状态协方差矩阵
    P01=(P1_IMM+(x1_IMM-x01)*(x1_IMM-x01)')*ui1(1)+...
        (P2_IMM+(x2_IMM-x01)*(x2_IMM-x01)')*ui1(2)+...
        (P3_IMM+(x3_IMM-x01)*(x3_IMM-x01)')*ui1(3);
    P02=(P1_IMM+(x1_IMM-x02)*(x1_IMM-x02)')*ui2(1)+...
        (P2_IMM+(x2_IMM-x02)*(x2_IMM-x02)')*ui2(2)+...
        (P3_IMM+(x3_IMM-x02)*(x3_IMM-x02)')*ui2(3);
    P03=(P1_IMM+(x1_IMM-x03)*(x1_IMM-x03)')*ui3(1)+...
        (P2_IMM+(x2_IMM-x03)*(x2_IMM-x03)')*ui3(2)+...
        (P3_IMM+(x3_IMM-x03)*(x3_IMM-x03)')*ui3(3); 
        
    %IMM基本算法
    [x1_IMM,P1_IMM,r1_IMM,S1_IMM]=Kalman(x01,P01,z_noise(:,:,t+1), F1,G,Q,H,R);
    [x2_IMM,P2_IMM,r2_IMM,S2_IMM]=Kalman(x02,P02,z_noise(:,:,t+1),F2,G,Q,H,R);
    [x3_IMM,P3_IMM,r3_IMM,S3_IMM]=Kalman(x03,P03,z_noise(:,:,t+1),F3,G,Q,H,R);
    %第三步--模型概率更新
    [u_IMM(:,t+1)]=Model_P_up(r1_IMM,r2_IMM,r3_IMM,S1_IMM,S2_IMM,S3_IMM,c_j);
    %第四步--模型综合
    [x_pro_IMM(:,t+1),P_IMM(:,:,t+1)]=Model_mix(x1_IMM,x2_IMM,x3_IMM,P1_IMM,P2_IMM,P3_IMM,u_IMM(:,t));


    %IMMPDA算法
    c_j=pij'*u_IMM_PDA(:,t);
    %计算每个模型的条件概率
    ui1=(1/c_j(1))*pij(:,1).*u_IMM_PDA(:,t);
    ui2=(1/c_j(2))*pij(:,2).*u_IMM_PDA(:,t);
    ui3=(1/c_j(3))*pij(:,3).*u_IMM_PDA(:,t);    
    
    % 计算各模型滤波初始化条件
    x01=x1_IMM_PDA*ui1(1)+x2_IMM_PDA*ui1(2)+x3_IMM_PDA*ui1(3);
    x02=x1_IMM_PDA*ui2(1)+x2_IMM_PDA*ui2(2)+x3_IMM_PDA*ui2(3);
    x03=x1_IMM_PDA*ui3(1)+x2_IMM_PDA*ui3(2)+x3_IMM_PDA*ui3(3);   
    
    %各模型滤波初始状态协方差矩阵
    P01=(P1_IMM_PDA+(x1_IMM_PDA-x01)*(x1_IMM_PDA-x01)')*ui1(1)+...
        (P2_IMM_PDA+(x2_IMM_PDA-x01)*(x2_IMM_PDA-x01)')*ui1(2)+...
        (P3_IMM_PDA+(x3_IMM_PDA-x01)*(x3_IMM_PDA-x01)')*ui1(3);
    P02=(P1_IMM_PDA+(x1_IMM_PDA-x02)*(x1_IMM_PDA-x02)')*ui2(1)+...
        (P2_IMM_PDA+(x2_IMM_PDA-x02)*(x2_IMM_PDA-x02)')*ui2(2)+...
        (P3_IMM_PDA+(x3_IMM_PDA-x02)*(x3_IMM_PDA-x02)')*ui2(3);
    P03=(P1_IMM_PDA+(x1_IMM_PDA-x03)*(x1_IMM_PDA-x03)')*ui3(1)+...
        (P2_IMM_PDA+(x2_IMM_PDA-x03)*(x2_IMM_PDA-x03)')*ui3(2)+...
        (P3_IMM_PDA+(x3_IMM_PDA-x03)*(x3_IMM_PDA-x03)')*ui3(3); 
        
    [x1_IMM_PDA,P1_IMM_PDA,L1] = PDA_Function(x01, P01, F1, H, z_noise(:,:,t+1));
    [x2_IMM_PDA,P2_IMM_PDA,L2] = PDA_Function(x02, P02, F2, H, z_noise(:,:,t+1));
    [x3_IMM_PDA,P3_IMM_PDA,L3] = PDA_Function(x03, P03, F3, H, z_noise(:,:,t+1));
    %第三步--模型概率更新
    [u_IMM_PDA(:,t+1)] = Model_P_up_PDA(L1,L2,L3,c_j);
    %第四步--模型综合
    [x_pro_IMM_PDA(:,t+1),P_IMM(:,:,t+1)]=Model_mix(x1_IMM_PDA,x2_IMM_PDA,x3_IMM_PDA,P1_IMM_PDA,P2_IMM_PDA,P3_IMM_PDA,u_IMM_PDA(:,t));   

    
    
%     [x4,P4,~,~] = Kalman(x4,P4,z_noise(:,:,t+1),F1,G,Q,H,R);
%     x4_rec(:,t+1) = x4;
end


figure
plot(x(1,:),x(3,:),'*-','LineWidth',1); hold on;grid on;
plot(x_pro_IMM(1,:),x_pro_IMM(3,:),'s-','LineWidth',1);
plot(x_pro_IMM_PDA(1,:),x_pro_IMM_PDA(3,:),'s-','LineWidth',1);
plot(x4_rec(1,:),x4_rec(3,:),'s-','LineWidth',1);
legend('Real','IMM','IMMPDA','CT');
% % 模型概率
% t=1:simTime;
% figure
% subplot(121)
% plot(t,u_IMM(1,t),'k.-',t,u_IMM(2,t),'r.-',t,u_IMM(3,t),'b.-');grid on
% title('IMM算法模型概率曲线');
% axis([0 simTime-1 0 1])
% xlabel('t/s'); ylabel('模型概率');
% legend('模型1','模型2','模型3');
% subplot(122)
% plot(t,u_IMM_PDA(1,t),'k.-',t,u_IMM_PDA(2,t),'r.-',t,u_IMM_PDA(3,t),'b.-');grid on
% title('IMM算法模型概率曲线');
% axis([0 simTime-1 0 1])
% xlabel('t/s'); ylabel('模型概率');
% legend('模型1','模型2','模型3');

% 位置误差
figure
subplot(2,1,1);
t=1:simTime;
plot(t,abs(x_pro_IMM(1,t)-x(1,t)),'LineWidth',1);hold on;grid on
plot(t,abs(x_pro_IMM_PDA(1,t)-x(1,t)),'LineWidth',1);hold on;grid on
% plot(t,abs(x4_rec(1,t)-x(1,t)),'LineWidth',1);grid on;
title('x坐标位置跟踪误差');
xlabel('t/s'); ylabel('x-error/m');
legend('IMM','IMMPDA','CT');

subplot(2,1,2);
t=1:simTime;
plot(t,abs(x_pro_IMM(3,t)-x(3,t)),'LineWidth',1);hold on;grid on
plot(t,abs(x_pro_IMM_PDA(3,t)-x(3,t)),'LineWidth',1);hold on;grid on
% plot(t,abs(x4_rec(3,t)-x(3,t)),'LineWidth',1);grid on;
title('y坐标位置跟踪误差');
xlabel('t/s'); ylabel('y-error/m');
legend('IMM','IMMPDA','CT');


%二、模型概率更新函数
function [u]=Model_P_up(r1,r2,r3,S1,S2,S3,c_j)
%模型概率更新函数
%计算似然函数
Lfun1=(1/sqrt(abs(2*pi*(det(S1)))))*exp((-1/2)*(r1'*inv(S1)*r1));   %Lfun1=1/(r1'*inv(S1)*r1);
Lfun2=(1/sqrt(abs(2*pi*(det(S2)))))*exp((-1/2)*(r2'*inv(S2)*r2));   %Lfun2=1/(r2'*inv(S2)*r2);
Lfun3=(1/sqrt(abs(2*pi*(det(S3)))))*exp((-1/2)*(r3'*inv(S3)*r3));    %Lfun3=1/(r3'*inv(S3)*r3);
c=[Lfun1,Lfun2,Lfun3]'.*c_j;
% 再归一化
u=[Lfun1,Lfun2,Lfun3]'.*c_j*(1/sum(c));

end

%二、模型概率更新函数
function [u]=Model_P_up_PDA(Lfun1,Lfun2,Lfun3,c_j)
%模型概率更新函数
c=[Lfun1,Lfun2,Lfun3]'.*c_j;
% 再归一化
u=[Lfun1,Lfun2,Lfun3]'.*c_j*(1/sum(c));

end

%三、状态混合函数
function [x_pro,P]=Model_mix(x1,x2,x3,P1,P2,P3,u)
x_pro=x1*u(1)+x2*u(2)+x3*u(3);
P=(P1+(x1-x_pro)*(x1-x_pro)')*u(1)+...
      (P2+(x2-x_pro)*(x2-x_pro)')*u(2)+...
      (P3+(x3-x_pro)*(x3-x_pro)')*u(3);
end

%四、Kalman滤波函数
function [X,P,e,S]=Kalman(X_Forward,P_Forward,z_noise,F,G,Q,H,R)
%卡尔曼滤波2012.2.27   IMM专用，参数略有不同
    x_predic = F*X_Forward;
    P_predic = F*P_Forward*F'+G*Q*G'; 
    Z_predic = H*x_predic;
    S = H*P_predic*H'+ R; 
    K = P_predic*H'*inv(S);
    y=[];m=0;
    for j=1:size(z_noise,1)
        d=z_noise(j,:)'-Z_predic;  
        D(j)=d'*inv(S)*d;  
    end
    [~,index] = min(D);
    Z = z_noise(index,:)';
    e = Z-Z_predic; %残差

    % 修正滤波值和误差协方差阵
    X=x_predic+K*e;
    M=K*H;
    n=size(M);
    I=eye(n);
    P=(I-K*H)*P_predic*(I-K*H)'+ K*R*K';
end


