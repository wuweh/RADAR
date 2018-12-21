% 参数设置
N = 100;   %粒子总数
Q = 5;      %过程噪声
R = 5;      %测量噪声
T = 20;     %测量时间
theta = pi/T;       %旋转角度
distance = 80/T;    %每次走的距离
WorldSize = 100;    %世界大小
X = zeros(2, T);    %存储系统状态
Z = zeros(2, T);    %存储系统的观测状态
P = zeros(2, N);    %建立粒子群
PCenter = zeros(2, T);  %所有粒子的中心位置
w = zeros(N, 1);         %每个粒子的权重
err = zeros(1,T);     %误差
X(:, 1) = [50; 20];     %初始系统状态
%wgn（m,n,p）产生一个m行n列的高斯白噪声的矩阵，p以dBW为单位指定输出噪声的强度。
Z(:, 1) = [50; 20] + wgn(2, 1, 10*log10(R));    %初始系统的观测状态

%初始化粒子群
for i = 1 : N
    P(:, i) = [WorldSize*rand; WorldSize*rand];%在worldSize区域内随机生成N个粒子
    dist = norm(P(:, i)-Z(:, 1));     %与测量位置相差的距离，用于估算该粒子的权重
    %由于上面已经随机生成了N个粒子，现在将其与真实的测量值z进行比较，越接近则权重越大，或者说差值越小权重越大
    %这里的权重计算是关于p(z/x)的分布，即观测方程的分布，假设观测噪声满足高斯分布，那么w（i）=p(z/x)
    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重
end
PCenter(:, 1) = sum(P, 2) / N;      %所有粒子的几何中心位置

%%
err(1) = norm(X(:, 1) - PCenter(:, 1));     %粒子几何中心与系统真实状态的误差
figure(1);
%利用set(gca,'propertyname','propertyvalue'......)命令可以调整图形的坐标属性。
%propertyname和property value分别为坐标属性名称（不区分字母大小写）及对应的属性值。
set(gca,'FontSize',10);
hold on%使当前及图形保持而不被刷新，准备接受此后将绘制的新曲线
plot(X(1, 1), X(2, 1), 'r.', 'markersize',30)   %系统状态位置
%axis[xmin,xmax,ymin,ymax]设定坐标范围，必须满足xmin<xmax，后面同理
axis([0 100 0 100]);
plot(P(1, :), P(2, :), 'kx', 'markersize',5);   %绘出各个粒子位置
plot(PCenter(1, 1), PCenter(2, 1), 'b.', 'markersize',25); %所有粒子的几何中心位置
%图形标注函数legend(string1,string2,string3...),在当前图中添加图例

legend('True State', 'Particles', 'The Center of Particles');
%添加图形标题命令title('string'),在当前坐标系的顶部加一个文本串string,作为该图形的标题
title('Initial State');
%使当前轴及图形不再具备不被刷新的性质
hold off

%%
%开始运动
for k = 2 : T

    %模拟一个弧线运动的状态
    %wgn（m,n,p）产生一个m行n列的高斯白噪声的矩阵，p以dBW为单位指定输出噪声的强度。
    X(:, k) = X(:, k-1) + distance * [(-cos(k * theta)); sin(k * theta)] + wgn(2, 1, 10*log10(Q));     %状态方程
    Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %观测方程 

    %粒子滤波
    %预测
    for i = 1 : N
        %将之前生成的粒子带入到状态方程，进行下一步状态的预测
        P(:, i) = P(:, i) + distance * [-cos(k * theta); sin(k * theta)] + wgn(2, 1, 10*log10(Q));
        dist = norm(P(:, i)-Z(:, k));     %与测量位置相差的距离，用于估算该粒子的权重
        %由于上面已经随机生成了N个粒子，现在将其与真实的测量值z进行比较，越接近则权重越大，或者说差值越小权重越大
        %这里的权重计算是关于p(z/x)的分布，即观测方程的分布，假设观测噪声满足高斯分布，那么w（i）=p(z/x)
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重
    end
%归一化权重
    wsum = sum(w);
    for i = 1 : N
        w(i) = w(i) / wsum;
    end

    %重采样（更新）
    for i = 1 : N
        wmax = 2 * max(w) * rand;  %另一种重采样规则,获取权重当中的最大值与均匀分布在0-1之间数的乘积再乘以2作为后面判断选出大权重粒子的依据
        %randi()函数生成均匀分布的伪随机整数，范围为imin--imax，如果没指定imin，则默认为1
        %r = randi([imin,imax],...)返回一个在[imin,imax]范围内的伪随机整数
        index = randi(N, 1);
        %在这里随机产生N个粒子当中的某个粒子，然后选择该粒子的权值，判断是否是大的粒子权重，是的话就把该粒子选出来
        while(wmax > w(index))
            %使vmax递减，不然的话单个粒子的权重是不可能大于vmax的值的
            wmax = wmax - w(index);
            index = index + 1;
            %如果从某个随机的粒子开始找，没找到从该粒子开始到最后的粒子权值总和大于vmax，那么就重新从第一个粒子开始查找
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     %从上面的index中得到新粒子
    end
    %sum(X,2)表示把X按行求和;如果是sum(X),那就是按列求和
    PCenter(:, k) = sum(P, 2) / N;      %所有粒子的几何中心位置

    %计算误差
    err(k) = norm(X(:, k) - PCenter(:, k));     %粒子几何中心与系统真实状态的误差

    figure(2);
    set(gca,'FontSize',12);
    %clf; 用来清除图形的命令。一般在画图之前用
    clf;
    hold on
    plot(X(1, k), X(2, k), 'r.', 'markersize',50);  %系统状态位置
    axis([0 100 0 100]);
    plot(P(1, :), P(2, :), 'k.', 'markersize',5);   %各个粒子位置
    plot(PCenter(1, k), PCenter(2, k), 'b.', 'markersize',25); %所有粒子的中心位置
    legend('True State', 'Particle', 'The Center of Particles');
    hold off
    pause(0.1);
end

%%
figure(3);
set(gca,'FontSize',12);
plot(X(1,:), X(2,:), 'r', Z(1,:), Z(2,:), 'g', PCenter(1,:), PCenter(2,:), 'b-');
axis([0 100 0 100]);
legend('True State', 'Measurement', 'Particle Filter');
xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);

%%
figure(4);
set(gca,'FontSize',12);
plot(err,'.-');
xlabel('t', 'FontSize', 20);
title('The err');