%使用说明：使用鼠标左键可以产生轨迹的各个点，鼠标右键为结束点。
axis([0 100 0 100])
hold on
% 程序变量初始化
xy = [];
n = 0;
% 使用循环，得到鼠标点击左键是的坐标位置
disp('Left mouse button picks points.')       %%在主窗口上提示操作方法
disp('Right mouse button picks last point.')
but = 1;
while but == 1
    [xi,yi,but] = ginput(1);         %%该函数可以得到鼠标左键的坐标点
    plot(xi,yi,'ro')
    n = n+1;
    xy(:,n) = [xi;yi];
end
% 利用插值函数获得光滑曲线，模拟目标运动曲线和测量数据。
t = 1:n;
ts = 1: 0.1: n;
xys = spline(t,xy,ts);
plot(xys(1,:),xys(2,:),'b-'); %%用蓝色直线画出目标运动曲线
Rx=10;Ry=10;    %%设置横轴、纵轴的测量噪声方差Rx、Ry。
plot(xys(1,:)+randn(size(xys(1,:)))*sqrt(Rx),xys(2,:)+randn(size(xys(1,:)))*sqrt(Ry),'k.');   %%，在2维坐标下用黑色点画出每一个测量数据
xlabel('x轴');ylabel('y轴')
hold off
%在另一张图上分别画出横、纵轴的目标运动数据和测量数据。
figure
subplot(2,1,1),plot(xys(1,:))
hold on
subplot(2,1,1),plot(xys(1,:)+randn(size(xys(1,:)))*sqrt(Rx),'.')
ylabel('x轴');
 
subplot(2,1,2),plot(xys(2,:))
hold on
subplot(2,1,2),plot(xys(2,:)+randn(size(xys(1,:)))*sqrt(Rx),'.')
ylabel('y轴')
hold off
%存储数据
save  mytarget1 xys ts