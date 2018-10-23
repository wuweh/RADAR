%三维Jerk模型------
%{
    
%}

clear all;
clc;
%定义参数,初始化
T=0.02;         %周期s
t0=T;           %初始时刻
t_end=10;    %终止时刻s
t=t0:T:t_end;      %采样时间段

N=length(t);  %采样时间点

x=zeros(1,N);     %x轴坐标初始化
y=zeros(1,N);     %y轴坐标初始化
z=zeros(1,N);     %z轴坐标初始化
x1=zeros(1,N);     %x轴速度初始化
y1=zeros(1,N);     %y轴速度初始化
z1=zeros(1,N);     %z轴速度初始化
x2=zeros(1,N);     %x轴加速度初始化
y2=zeros(1,N);     %y轴加速度初始化
z2=zeros(1,N);     %z轴加速度初始化
x3=zeros(1,N);     %x轴加加速度初始化
y3=zeros(1,N);     %y轴加加速度初始化
z3=zeros(1,N);     %z轴加加速度初始化

xdelta_R=20;       %雷达x方向测距误差m
ydelta_R=30;       %雷达y方向测距误差m
zdelta_R=25;       %雷达z方向测距误差m

alpha=0.1;
jerk_x=-20;
g=-9.8;
X=[x;  x1; x2;  x3;    y;  y1;  y2;  y3;    z;  z1;  z2; z3];       %目标状态矩阵
%赋初值
X(:,1)=[3000;  -200;  0;  jerk_x;    4000;  200;  0;  0;    10000;  0;  g; 0];
%{
x(1)=3000;          %x初始位置m
y(1)=4000;          %y初始位置m
vx(1)=200;          %x初始速度m/s
vy(1)=200;          %y初始速度m/s
vz(1)=0;            %z初始速度m/s
%}

p1=(2 - 2*alpha*T + alpha^2*T^2 - 2*exp(-alpha*T))/(2*alpha^3);
q1=(exp(-alpha*T) - 1 + alpha*T)/alpha^2;
r1=(1 - exp(-alpha*T))/alpha;
s1=exp(-alpha*T);

%目标的状态更新
F4=[ 1, T, T^2/2, p1;
     0, 1, T,     q1;
     0, 0, 1,     r1;
     0, 0, 0,     s1 ];
 
F=[F4,              zeros(4,4),     zeros(4,4);
   zeros(4,4),      F4,             zeros(4,4);
   zeros(4,4),      zeros(4,4),     F4];

 H=[ 1, 0, 0, 0,  0, 0, 0,  0,  0, 0, 0, 0;
        0, 0, 0,  0, 1, 0, 0, 0,   0, 0, 0, 0;
        0, 0, 0, 0,  0, 0, 0, 0,   1, 0, 0, 0  ];
Y=zeros(size(H,1),N);   %初始化测量状态.
r=0.05;
for k=2:N 
    X(:,k)=F*X(:,k-1);
    X(11,k)=g;
    Y(:,k)=H*X(:,k-1)+r*randn(size(H,1),1);
end
    
x=X(1,:);
y=X(5,:);
z=X(9,:);

for i=1:1:N
    R(i)=sqrt( X(1,i)^2 + X(5,i)^2 + X(9,i)^2 );
    v(i)=sqrt( X(2,i)^2 + X(6,i)^2 + X(10,i)^2 );
    a(i)=sqrt( X(3,i)^2 + X(7,i)^2 + X(11,i)^2 );
end

figure(1)
plot3(x,y,z);
xlabel('x轴:m'); zlabel('y轴:m');grid on;  

figure(2)
subplot(311),plot(t,R,'-b');
legend('R'),title('R曲线图'); xlabel('时间：s'); ylabel('R:m');
subplot(312),plot(t,v);
legend('v'),title('速度曲线'); xlabel('采样时间点'); ylabel('速度:m/s');
subplot(313),plot(t,a);
legend('a'),title('加速度曲线'); xlabel('t:s'); ylabel('加速度:m/s^2');
