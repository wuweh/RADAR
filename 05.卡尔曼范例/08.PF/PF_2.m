%%自由落体过程的粒子滤波
clear all;close all;clc
% 参数设置
N = 200;   %粒子总数
t = 20;                     % 仿真时间
Ts = 0.1;                   % 采样周期 
len = fix(t/Ts);            % 仿真步数
kx = .01;   ky = .05;       % 阻尼系数
g = 9.8;                    % 重力
X(1,:) = [0, 50, 500, 0]; % 状态模拟的初值
dax = 2; day = 2;       % 系统噪声

%生成真实轨迹
for k=2:len
    x = X(k-1,1); vx = X(k-1,2); y = X(k-1,3); vy = X(k-1,4); 
    x = x + vx*Ts;
    vx = vx + (-kx*vx^2+dax*randn)*Ts;
    y = y + vy*Ts;
    vy = vy + (ky*vy^2-g+day*randn)*Ts;
    X(k,:) = [x, vx, y, vy];
end

%生成量测
mrad = 0.001;
dr = 4; 
dafa = 10*mrad; % 量测噪声
for k=1:len
    r = sqrt(X(k,1)^2+X(k,3)^2) + dr*randn(1,1);
    a = atan(X(k,1)/X(k,3)) + dafa*randn(1,1);
    Z(k,:) = [r, a];
end

V = [1,0.01,1,0.01];
% 用一个高斯分布随机的产生初始的粒子  
for i = 1:N    
    x_P(i,:) = X(1,:) + sqrt(V) * randn;    
end

PCenter(1,:) = mean(x_P); 

%开始运动
for k = 2: len
    %粒子滤波
    %预测
    %对每一个粒子进行一次状态更新
    for i = 1 : N
        x = x_P(i,1); vx = x_P(i,2); y = x_P(i,3); vy = x_P(i,4);
        x = x + vx*Ts;
        vx = vx + (-kx*vx^2+dax*randn(1,1))*Ts;
        y = y + vy*Ts;
        vy = vy + (ky*vy^2-g+day*randn(1))*Ts;
        x_P_update(i,:) = [x, vx, y, vy];
        dist = sqrt((x_P_update(i,1)-Z(k,1)*sin(Z(k,2)))^2+(x_P_update(i,3)-Z(k,1)*cos(Z(k,2)))^2);       
        w(i) = (1/sqrt(dr)/sqrt(2*pi)) * exp(-(dist)^2/2/dr);   %求权重
    end
    
    %归一化权重
    w = w./sum(w);

    %重采样（方法一）
%     for i = 1 : N
%         wmax = 2 * max(w) * rand;  %另一种重采样规则,获取权重当中的最大值与均匀分布在0-1之间数的乘积再乘以2作为后面判断选出大权重粒子的依据
%         %randi()函数生成均匀分布的伪随机整数，范围为imin--imax，如果没指定imin，则默认为1
%         %r = randi([imin,imax],...)返回一个在[imin,imax]范围内的伪随机整数
%         index = randi(N, 1);
%         %在这里随机产生N个粒子当中的某个粒子，然后选择该粒子的权值，判断是否是大的粒子权重，是的话就把该粒子选出来
%         while(wmax > w(index))
%             %使vmax递减，不然的话单个粒子的权重是不可能大于vmax的值的
%             wmax = wmax - w(index);
%             index = index + 1;
%             %如果从某个随机的粒子开始找，没找到从该粒子开始到最后的粒子权值总和大于vmax，那么就重新从第一个粒子开始查找
%             if index > N
%                 index = 1;
%             end          
%         end
%         x_P(i,:) = x_P_update(index,:);     %从上面的index中得到新粒子
%     end
    
    %重采样（方法二）
    for i = 1 : N    
        x_P(i,:) = x_P_update(find(rand <= cumsum(w),1),:);   % 粒子权重大的将多得到后代    
    end  

    PCenter(k,:) = mean(x_P); 
    %计算误差
    err(k) = sqrt((PCenter(k,1)-X(k,1))^2 + (PCenter(k,3)-X(k,3))^2);     %粒子几何中心与系统真实状态的误差

%     figure(1);
%     set(gca,'FontSize',12);
    plot(X(k,1), X(k,3), 'r*', 'markersize',3);    hold on  %系统状态位置
    plot(x_P(:,1), x_P(:,3), 'k.', 'markersize',5);     hold on  %各个粒子位置
    plot(PCenter(k,1), PCenter(k,3), 'b*', 'markersize',3);    hold on ;grid on%所有粒子的中心位置
    axis([0 300 0 500])
%     legend('True State', 'Particle', 'The Center of Particles');
     pause(0.001);
end

%作图部分
figure;
set(gca,'FontSize',12);
plot(X(:,1), X(:,3), 'b.-');hold on;
plot(Z(:,1).*sin(Z(:,2)), Z(:,1).*cos(Z(:,2)),'.');hold on;
plot(PCenter(:,1), PCenter(:,3), 'r.-');hold on;
grid on;axis([0 350 0 550])
legend('True State', 'Measurement', 'Particle Filter');
xlabel('x', 'FontSize', 10); ylabel('y', 'FontSize', 10);

figure;
set(gca,'FontSize',10);
plot(err,'.-');
xlabel('t', 'FontSize', 10);
title('The err');