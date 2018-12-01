%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  给定α、β系数条件下，验证常增益α-β滤波（此系数与观测信号噪声有关，具体参考文献）
%  输入：真实数据（位置+速度）realval
%            观测数据（位置）obseval 加观测噪声
%  输出：估计数据（位置+速度）estimval
%            校正输出值（位置+速度）reviseval
%  参考文献：蔡庆宇、张伯彦，曲洪泉.相控阵雷达数据处理教程.北京：电子工业出版社.pp.60-66
%  樊文贵/电子信息工程学院/北京航空航天大学
%  5/8/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;
clear all;

alpha=0.27;
beta=0.04229;

T=0.1;                      %采样间隔/s
phi=[1,T;0,1];              %转移矩阵
K=[alpha;beta/T];           %常数增益矩阵
H=[1,1];                    %测量矩阵

%% 轨道生成(输入观测值obseval)
t=0:T:100;                  %仿真时间/s
LEN=length(t);              %仿真数据长度

xi=0;                       %位置观测噪声均值
zeta=40e3;                  %位置观测噪声方差
realval=zeros(2,LEN);                                 

realval(1,:)=0.3e3*t+100e3;         %位置真实值/m
realval(2,:)=0.3e3*ones(1,LEN);     %速度真实值/m/s

origsite=100e3;                     %初始位置/m
origvelo=0.3e3;                     %初始速度/m/s

obseval=H*realval+(sqrt(zeta)*randn(1,LEN)+xi*ones(1,LEN)); %位置观测值 +方差为zeta均值为0的高斯噪声

figure(1);
subplot(221);plot(t,realval(1,:)/1e3);title('目标真实轨道');
xlabel('时间/s');ylabel('距离/km');grid on;

%% 轨道预测与校正
estimval=zeros(2,LEN);       %轨道预测值(估计值)
reviseval=zeros(2,LEN);      %轨道校正值(输出值)

estimval(1,1)=origsite;        %初始化估计值（位置/m）
estimval(2,1)=origvelo;       %初始化估计值（速度/m/s）

reviseval(1,1)=origsite;        %初始化校正值（位置/m）
reviseval(2,1)=origvelo;       %初始化校正值（速度/m/s）

%%分析对象有两个：一个是距离；一个是速度
for n=1:1:LEN-1
    %由k时刻校正输出值预测k+1时刻状态
    %estimval是卡尔曼模型中的预测值
    %phi是转移矩阵
    estimval(:,n+1)=phi*reviseval(:,n);
    
    %在k+1时刻获得实际观测值obseval，在此基础上校正输出值
    %reviseval是卡尔曼模型中的最优输出值，其中K为卡尔曼增益系数，在α-β滤波模型中，为一个定值
    %H是量测矩阵
    reviseval(:,n+1)=estimval(:,n+1)+K*(obseval(:,n+1)-H*estimval(:,n+1));
end

%% 输出结果
subplot(222);
plot(t,estimval(1,:)/1e3);title('估计轨道');
xlabel('时间/s');ylabel('距离/km');grid on;

subplot(223);
plot(t,reviseval(1,1:LEN)/1e3);title('校正(输出)轨道');
xlabel('时间/s');ylabel('距离/km');grid on;

subplot(224);
plot(t,abs(estimval(1,1:LEN)-reviseval(1,1:LEN))/1e3);title('误差');
xlabel('时间/s');ylabel('误差/km');grid on;

%位置
figure;
plot(t,estimval(1,1:LEN)/1e3);hold on;
plot(t,reviseval(1,1:LEN)/1e3,'r');hold on;
plot(t,realval(1,1:LEN)/1e3,'k');hold off;legend('估计值','校正值','观测值');grid on;

%速度
figure;
plot(t,estimval(2,1:LEN));hold on;
plot(t,reviseval(2,1:LEN),'r');hold on;
plot(t,realval(2,1:LEN),'k');hold off;legend('估计值','校正值','观测值');