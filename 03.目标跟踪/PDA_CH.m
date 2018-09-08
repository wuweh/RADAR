% PDA算法仿真代码
% 在原始的PDA仿真前提下，加入了对X轴的平滑

clc
clear all
close all

tic
T = 1;
n=50;
global Pd Pg gamma C;
Pd=1;       %检测概率 
Pg=0.99;    %正确量测落入跟踪门内得概率 

g_sigma = 9.21;     %门限 
gamma = 1;          %每一个单位面积内产生一个杂波 ，虚假测量的空间密度

target_position=zeros(4,n);   %目标观测互联存储矩阵    

x_filter=zeros(4,n); 
x_filter1=zeros(4,n);

A = [1 T 0 0;
     0 1 0 0;
     0 0 1 T;
     0 0 0 1];  %状态转移矩阵--->假设做匀速直线运动
 
Q = [ 0.01 0 0 0; 
      0 0.01 0 0; 
      0 0 0.01 0; 
      0 0 0 0.01];     %实际过程噪声协方差矩阵
 
C = [1 0 0 0;
     0 1 0 0
     0 0 1 0;
     0 0 0 1];     %观测矩阵：得到 X Y的量测值

target_delta=1;           %目标对应的观测标准差 
R=[ 1^2 0 0 0;
    0 0.1^2 0 0;
    0 0 1^2 0 ;
    0 0 0 0.1^2]; %观测协方差

% G=[T^2/2 0;
%    0 T^2/2;
%    T 0;
%    0 T]; 

R11=target_delta; 
R22=target_delta; 
R12=0; 
R21=0;

P=[R11 R11/T R12 R12/T;     R11/T 2*R11/T^2 R12/T 2*R12/T^2; 
   R21 R21/T R22 R22/T;     R21/T 2*R21/T^2 R22/T 2*R22/T^2];   %初始协方差 

%初始状态
X0=[200;0;500;15];            
target_position(:,1)=X0; 
Vk=[target_delta*randn;
    0.1^2;
    target_delta*randn;
    0.1^2]; 
Zk(:,1)=C*target_position(:,1)+Vk;   

%************************************************
%          量测生成
%************************************************
for i=2:1:n 
    target_position(:,i)=A*target_position(:,i-1); % 循环获得50个真实状态
    
%     Vk=[target_delta*randn;target_delta*randn];   %生成量测噪声
    Vk=[target_delta*randn;
        0.1^2;
        target_delta*randn;
        0.1^2]; 
    Zk(:,i)=C*target_position(:,i)+Vk;      %生成50个量测值 
end

Noise=[];
NOISE=[]; 

sumx = 0;
figure;
%滤波开始 
for t=1:n
    if t~=1 
        x_predic = A*x_filter(:,t-1);   %用前一时刻的滤波值来预测当前的值  
    else 
        x_predic = target_position(:,1); %第一次采样我们用真实位置当预测值  
    end 
    
    plot(Zk(1,t),Zk(3,t),'g*') ;
    hold on; 

    %标准卡尔曼滤波过程
    %详见《目标跟踪前沿理论与应用》P42
    %计算预测协方差和预测值
    P_predic = A*P*A'+Q; 
    Z_predic = C*x_predic;
    
    %计算量测协方差（新息协方差）
    S = C*P_predic*C'+ R; 
    K = P_predic*C'*inv(S); %增益 
    
    %% 下列代码模拟实际杂波环境
    %雷达数据处理及应用 P52 (8.53)
    Av=pi*g_sigma*sqrt(det(S))   ;                         %计算椭球体积，这里算的是面积    
%      number_returns=floor(10*Av*gamma+1) ;                 %错误回波数 
     number_returns = 20;
    %雷达数据处理及应用 P127 (7.49)
    side=sqrt(10*Av)/2;                                     %求出正方行边长的二分之一 
    %雷达数据处理及应用 P127 (7.46)
    Noise_x=x_predic(1)-side+0.1*rand(1,number_returns)*side;            %在预测值周围产生多余回波 
    Noise_y=x_predic(3)-side+0.1*rand(1,number_returns)*side; 
    Noise_vx= x_predic(2)-side+0.1*rand(1,number_returns)*side;            %在预测值周围产生多余回波 
    Noise_vy= x_predic(4)-side+0.1*rand(1,number_returns)*side; 
    Noise=[Noise_x ;Noise_vx;Noise_y;Noise_vy]; 
     
    b=zeros(1,4); 
    b(1)=Zk(1,t); 
    b(2)=Zk(2,t); 
    b(3)=Zk(3,t); 
    b(4)=Zk(4,t); 
    y1=[Noise b'];   %将接收到的所有的回波存在y1中 
    y=[]; 
    d=[]; 
    m=0; 
    
    %% 计算落入门限内的杂波，此处采样了卡方分布算法
    for j=1:(number_returns+1)     
        d=y1(:,j)-Z_predic; 
        D=d'*inv(S)*d;  %卡方分布值计算
        if D<=g_sigma 
            plot(y1(1,j),y1(3,j),'.');
            hold on;
            y=[y y1(:,j)];   %把落入跟踪门中的所有回波放入y中 
            m=m+1;          %记录观测个数 
        end 
    end 
    
    %% PDA_Function
    [x_putput, P] = PDA_Function(x_predic, P_predic, S, Z_predic, K,y, m);
    
    x_filter(:,t) = x_putput;
    sumx = sumx + x_filter(1,t);
    x_filter1(:,t) = x_filter(:,t);
    x_filter1(1,t) = sumx/t;
    
end

plot(target_position(1,:),target_position(3,:),'-bo');
hold on; 
plot(x_filter(1,:),x_filter(3,:),'-r*');
hold on; 
plot(x_filter1(1,:),x_filter1(3,:),'-ro');
hold off; 
% legend('观测位置','杂波位置','真实位置','预测位置','平滑位置');
axis([190 210 400 1400]);
toc

