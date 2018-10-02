%% ----------------------------------------------------------------
%�ݻ��������˲���CKF��������������ƽ������ģ��
%
%% ----------------------------------------------------------------
function CubatureKF
clear all;
close all;
clc;


n=1;%ϵͳ��ά��
m=2*n;%�ݻ�����
w=1/m;%Ȩֵw=1/m
kesi=sqrt(m/2)*[1,-1];%kesi=sqrt(m/2)*[1]_i
Q=10;%��������
R=1;%��������

x=0.1;
Pplus=10;
xhat=x;%x(0|0)�ĳ�ʼ-ֵ���ȡֵ
xarray=[x];
zarray=[x^2/20+sqrt(R)*randn];
xhatarray=[x];

num=100;%���泤��
for i=1:num
    x=0.5*x+25*x/(1+x^2)+8*cos(1.2*(i-1))+sqrt(Q)*randn;%ϵͳ����
    z=x^2/20+sqrt(R)*randn;%���ⷽ��
    xarray=[xarray x];
    zarray=[zarray,z];
%% ----------------------------CKF�˲�----------------------------

%% ----------------------------ʱ�����----------------------------
    %��1��Э�������Cholesky�ֽ�
    Shat=chol(Pplus,'lower');
    for cpoint=1:m
        %��2�������ݻ���
        rjpoint(cpoint)=Shat*kesi(cpoint)+xhat;
        %��3�������ݻ���
        Xminus(cpoint)=0.5*rjpoint(cpoint)+25*rjpoint(cpoint)/(1+rjpoint(cpoint)^2)+8*cos(1.2*(i-1)); %�ݻ��㾭�������Ժ������ֵ
    end
    %��4��״̬Ԥ��
    xhat=w*sum(Xminus);
    %��5��״̬Ԥ��Э������
    Pminus=w*sum(Xminus.^2)-xhat*xhat'+Q;
%% ---------------------------------------------------------------

%% ----------------------------�������----------------------------
    %��1������Cholesky�ֽ�
    Sminus=chol(Pminus,'lower');
    for cpoint=1:m
        %��2�������ݻ���
        rjpoint1(cpoint)=Sminus*kesi(cpoint)+xhat;
        %��3�������ݻ���
        Z(cpoint)=rjpoint1(cpoint)^2/20;%�ݻ��㾭�������Ժ������ֵ
    end    
    %��4���۲�Ԥ��
    zhat=w*sum(Z);
    %��5���۲�Ԥ��Э������
    Pzminus=w*sum(Z.^2)-zhat^2+R;
    %��6����Э������
    Pxzminus=w*rjpoint1*Z'-xhat*zhat;
    %��7�����㿨��������
    W=Pxzminus*inv(Pzminus);
    %��8��״̬����
    xhat=xhat+W*(z-zhat);
    %��9��״̬Э����������
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
