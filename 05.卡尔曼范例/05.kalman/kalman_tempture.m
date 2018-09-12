clear
clc;
N=300;
CON = 25;%房间温度，假定温度是恒定的

%%%%%%%%%%%%%%%kalman filter%%%%%%%%%%%%%%%%%%%%%%
x = zeros(1,N);
y = 2^0.5 * randn(1,N) + CON;%加过程噪声的状态输出

x(1) = 1;
p = 10;

Q = cov(randn(1,N));%过程噪声协方差
R = cov(randn(1,N));%观测噪声协方差
for k = 2 : N
    x(k) = x(k - 1);    %预估计k时刻状态变量的值
    p = p + Q;  %预估值的协方差
    kg = p / (p + R);%kalman gain  卡尔曼增益
    
    %由预计值x(k)和测量值z(k)以及卡尔曼增益kg得到最优输出值
    x(k) = x(k) + kg * (y(k) - x(k));
   
    p = (1 - kg) * p; %更新最优估计值的协方差
end


%%%%%%%%%%%Smoothness Filter%%%%%%%%%%%%%%%%%%%%%%%%

Filter_Wid = 10;
smooth_res = zeros(1,N);
for i = Filter_Wid + 1 : N
tempsum = 0;
for j = i - Filter_Wid : i - 1
tempsum = tempsum + y(j);
end
smooth_res(i) = tempsum / Filter_Wid;
end
% figure(1);
% hist(y);
t=1:N;
figure(1);

expValue = zeros(1,N);
for i = 1: N
    expValue(i) = CON;
end
plot(t,expValue,'r',t,x,'g',t,y,'b',t,smooth_res,'k');
legend('expected','estimate','measure','smooth result');
axis([0 N 20 30])
xlabel('Sample time');
ylabel('Room Temperature');
title('Smooth filter VS kalman filter');