%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  功能说明： UKF在目标跟踪中的应用
%  参数说明： 1、状态6维，x方向的位置、速度、加速度；
%                y方向的位置、速度、加速度；
%             2、观测信息为距离和角度；
%  详细原理介绍及中文注释请参考：
%  《卡尔曼滤波原理及应用-MATLAB仿真》，电子工业出版社，黄小平著。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ukf_for_track_6_div_system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=6;
t=0.05;
%过程噪声协方差
Q=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 0.01 0 0 0;
    0 0 0 0.01 0 0;
    0 0 0 0 0.0001 0;
    0 0 0 0 0 0.0001];
%测量噪声协方差
R = [100 0;
    0 0.001^2];
%状态方程
%x1为x轴位置，x2为y轴位置，x3、x4分别是x、y轴的速度，x5、x6分别是加速度
f=@(x)[ x(1)+t*x(3)+0.5*t^2*x(5);
        x(2)+t*x(4)+0.5*t^2*x(6); 
        x(3)+t*x(5);
        x(4)+t*x(6);
        x(5);
        x(6)];
%观测方程
h=@(x)[ sqrt(x(1)^2+x(2)^2);
        atan(x(2)/x(1))];

s1=[92000;82000;-130;-230;0;0];%实际位置初始化
s2=[91000;82000;-130;-230;0;0];%实际位置初始化
s3=[90000;82000;-130;-230;0;0];%实际位置初始化
s4=[89000;82000;-130;-230;0;0];%实际位置初始化
s5=[90550;90930;-250;190;32;-23];%实际位置初始化
s6=[90120;90600;-139;260;47;-24];%实际位置初始化
s7=[90990;91390;-200;320;53;-25];%实际位置初始化



x0=s1+sqrtm(Q)*randn(n,1);%先验位置初始化状态

%初始化协方差
P0 =[100 0 0 0 0 0;
    0 100 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 0.0001 0;
    0 0 0 0 0 0.0001];

%总仿真步数
N=50;

%ukf滤波状态初始化
Xukf = zeros(n,N);
X1 = zeros(n,N);
X2 = zeros(n,N);
X3 = zeros(n,N);
X4 = zeros(n,N);
X5 = zeros(n,N);
X6 = zeros(n,N);
X7 = zeros(n,N);
Z = zeros(2,N);

figure
hold on;
for i=1:N
    X1(:,i)= f(s1)+sqrtm(Q)*randn(6,1); %模拟，产生目标运动真实轨迹
    s1 = X1(:,i);
     s1_1(i) = sqrt(X1(1,i)^2+X1(2,i)^2);
    
    X2(:,i)= f(s2)+sqrtm(Q)*randn(6,1); %模拟，产生目标运动真实轨迹
    s2 = X2(:,i);
     s2_1(i) = sqrt(X2(1,i)^2+X2(2,i)^2);
    
    X3(:,i)= f(s3)+sqrtm(Q)*randn(6,1); %模拟，产生目标运动真实轨迹
    s3 = X3(:,i);
     s3_1(i) = sqrt(X3(1,i)^2+X3(2,i)^2);
    
    X4(:,i)= f(s4)+sqrtm(Q)*randn(6,1); %模拟，产生目标运动真实轨迹
    s4 = X4(:,i);
     s4_1(i) = sqrt(X4(1,i)^2+X4(2,i)^2);
    
    X5(:,i)= f(s5)+sqrtm(Q)*randn(6,1); %模拟，产生目标运动真实轨迹
    s5 = X5(:,i);
     s5_1(i) = sqrt(X5(1,i)^2+X5(2,i)^2);
    
    X6(:,i)= f(s6)+sqrtm(Q)*randn(6,1); %模拟，产生目标运动真实轨迹
    s6 = X6(:,i);
     s6_1(i) = sqrt(X6(1,i)^2+X6(2,i)^2);
    
    X7(:,i)= f(s7)+sqrtm(Q)*randn(6,1); %模拟，产生目标运动真实轨迹
    s7 = X7(:,i);
     s7_1(i) = sqrt(X7(1,i)^2+X7(2,i)^2); 
     
    plot(i,s1_1(i),'.g');
    plot(i,s2_1(i),'.r');
    plot(i,s3_1(i),'.y');
    plot(i,s4_1(i),'.b');
    plot(i,s5_1(i),'*g');
    plot(i,s6_1(i),'*r');
    plot(i,s7_1(i),'*y');
end
figure;
hold on;
plot(X1(1,:),X1(2,:),'.g');
plot(X2(1,:),X2(2,:),'.r');
plot(X3(1,:),X3(2,:),'.y');
plot(X4(1,:),X4(2,:),'.b');
plot(X5(1,:),X5(2,:),'*g');
plot(X6(1,:),X6(2,:),'*r');
plot(X7(1,:),X7(2,:),'*y');

ux=x0;
for k=1:N
    %h(X(:,k))指的是通过测量得到目标的位置信息：距离和角度
    Z(:,k)= h(X1(:,k)) + sqrtm(R)*randn(2,1);        %测量值加上噪声
    [Xukf(:,k), P0] = ukf(f,ux,P0,h,Z(:,k),Q,R);    %调用ukf算法
    ux=Xukf(:,k);
end

%跟踪误差分析
for k=1:N
    RMS(k)=sqrt((X1(1,k)-Xukf(1,k))^2+(X1(2,k)-Xukf(2,k))^2 );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
t=1:N;
hold on;box on;
plot(X1(1,t),X1(2,t), 'k-')  %真实目标xy信息
plot(Z(1,t).*cos(Z(2,t)),Z(1,t).*sin(Z(2,t)),'-b.')%测量目标信息
plot(Xukf(1,t),Xukf(2,t),'-r.')%ukf目标信息
legend('实际值','测量值','ukf估计值');
xlabel('x方向位置/米');
ylabel('y方向位置/米');

figure
box on;
plot(RMS,'-ko','MarkerFace','r')
xlabel('t/秒')
ylabel('偏差/米')

%%
%ffun：状态转移矩阵
%X：先验状态
%P：协方差
% hfun:观测方程
% Z:当前观测状态
% Q:过程噪声协方差
% R:测量噪声协方差
function [X,P]=ukf(ffun,X,P,hfun,Z,Q,R)
L=numel(X);  %获得X的元素个数
m=numel(Z);  %获得Z的元素个数
%三个参数设置
alpha=1e-2;
ki=0;
beta=2;

lambda=alpha^2*(L+ki)-L; %公式1
c=L+lambda;              %中间变量
Wm=[lambda/c 0.5/c+zeros(1,2*L)]; %公式2
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);%公式3

c=sqrt(c);
Xsigmaset=sigmas(X,P,c); 
[X1means,X1,P1,X2]=ut(ffun,Xsigmaset,Wm,Wc,L,Q);   
[Zpre,Z1,Pzz,Z2]=ut(hfun,X1,Wm,Wc,m,R);

Pxz=X2*diag(Wc)*Z2';
K=Pxz/(Pzz); %滤波增益
X=X1means+K*(Z-Zpre);   %后验状态均值
P=P1-K*Pxz';            %后验状态协方差

%%
function [Xmeans,Xsigma_pre,P,Xdiv]=ut(fun,Xsigma,Wm,Wc,n,COV)
LL=size(Xsigma,2);
Xmeans=zeros(n,1);
Xsigma_pre=zeros(n,LL);
for k=1:LL
    Xsigma_pre(:,k)=fun(Xsigma(:,k));  %状态的sigma点非线性变换
    Xmeans=Xmeans+Wm(k)*Xsigma_pre(:,k); %状态均值
end
Xdiv=Xsigma_pre-Xmeans(:,ones(1,LL));
P=Xdiv*diag(Wc)*Xdiv'+COV;              %状态协方差

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sigma points around reference point
% 构造2n+1个sigma点
% Inputs:
% x: reference point
% P: covariance
% c: coefficient
% Output:
% X: Sigma points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xset=sigmas(X,P,c)
A = c*chol(P)';
Y = X(:,ones(1,numel(X)));
Xset = [X Y+A Y-A];
