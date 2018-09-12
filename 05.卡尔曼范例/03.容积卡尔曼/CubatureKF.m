%% ----------------------------------------------------------------
%容积卡尔曼滤波（CKF）——单变量非平稳增长模型
%
%% ----------------------------------------------------------------
function CubatureKF
clear all;
close all;
clc;


n=1;%系统的维数
m=2*n;%容积点数
w=1/m;%权值w=1/m
kesi=sqrt(m/2)*[1,-1];%kesi=sqrt(m/2)*[1]_i
Q=10;%过程噪声
R=1;%量测噪声

x=0.1;
Pplus=10;
xhat=x;%x(0|0)的初始-值随机取值
xarray=[x];
zarray=[x^2/20+sqrt(R)*randn];
xhatarray=[x];

num=100;%仿真长度
for i=1:num
    x=0.5*x+25*x/(1+x^2)+8*cos(1.2*(i-1))+sqrt(Q)*randn;%系统方程
    z=x^2/20+sqrt(R)*randn;%量测方程
    xarray=[xarray x];
    zarray=[zarray,z];
%% ----------------------------CKF滤波----------------------------

%% ----------------------------时间更新----------------------------
    %（1）协方差矩阵Cholesky分解
    Shat=chol(Pplus,'lower');
    for cpoint=1:m
        %（2）计算容积点
        rjpoint(cpoint)=Shat*kesi(cpoint)+xhat;
        %（3）传播容积点
        Xminus(cpoint)=0.5*rjpoint(cpoint)+25*rjpoint(cpoint)/(1+rjpoint(cpoint)^2)+8*cos(1.2*(i-1)); %容积点经过非线性函数后的值
    end
    %（4）状态预测
    xhat=w*sum(Xminus);
    %（5）状态预测协方差阵
    Pminus=w*sum(Xminus.^2)-xhat*xhat'+Q;
%% ---------------------------------------------------------------

%% ----------------------------量测更新----------------------------
    %（1）矩阵Cholesky分解
    Sminus=chol(Pminus,'lower');
    for cpoint=1:m
        %（2）计算容积点
        rjpoint1(cpoint)=Sminus*kesi(cpoint)+xhat;
        %（3）传播容积点
        Z(cpoint)=rjpoint1(cpoint)^2/20;%容积点经过非线性函数后的值
    end    
    %（4）观测预测
    zhat=w*sum(Z);
    %（5）观测预测协方差阵
    Pzminus=w*sum(Z.^2)-zhat^2+R;
    %（6）互协方差阵
    Pxzminus=w*rjpoint1*Z'-xhat*zhat;
    %（7）计算卡尔曼增益
    W=Pxzminus*inv(Pzminus);
    %（8）状态更新
    xhat=xhat+W*(z-zhat);
    %（9）状态协方差矩阵更新
    Pplus=Pminus-W*Pzminus*W';
%% ---------------------------------------------------------------

    xhatarray=[xhatarray xhat];    
end
%% ---------------------------------------------------------------


k=0:num;
figure(1)
plot(k,xarray,'b.',k,xhatarray,'r-');
set(gca,'fontname','Times New Roman','fontsize',12);
set(gcf,'Color','White');
xlabel('Time step','fontname','Times New Roman','fontsize',16);
ylabel('State','fontname','Times New Roman','fontsize',16);
axis tight;
legend('True state','CKF estimates');
title('CKF estimates','fontname','Times New Roman','fontsize',16) ;

error=xarray-xhatarray;
CKF_RMS=rms(error);
fa(num)=CKF_RMS;

figure(2)
plot(k,error);
set(gca,'fontname','Times New Roman','fontsize',12);
set(gcf,'Color','White');
xlabel('Time step','fontname','Times New Roman','fontsize',16); 
ylabel('State','fontname','Times New Roman','fontsize',16);
axis tight;
title('CKF estimate errors','fontname','Times New Roman','fontsize',16) ;

disp(['Cubature Kalman filter RMS error = ', num2str(CKF_RMS)]);
