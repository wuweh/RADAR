%��չKalman�˲���Ŀ������е�Ӧ��
%function EKF_For_TargetTracking
clc;clear;
T=1;%�״�ɨ������

N=60/T;%�ܵĲ�������

X=zeros(4,N);%Ŀ����ʵλ�á��ٶ�

X(:,1)=[-100,2,200,20];%Ŀ���ʼλ�á��ٶ�
 	
Z=zeros(1,N);%��������λ�õĹ۲�

delta_w=1e-3;%����������������Ŀ�����ʵ�켣����������

Q=delta_w*diag([0.5,1]);%������������

G=[T^2/2,0;T,0;0,T^2/2;0,T];%����������������
	
R=5;%�۲���������
 
F=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];%״̬ת�ƾ���

x0=200;%�۲�վ��λ��

y0=300;
 	
Xstation=[x0,y0];

for t=2:N
    X(:,t)=F*X(:,t-1)+G*sqrtm(Q)*randn(2,1);%Ŀ����ʵ�켣	
end

for t=1:N	
    Z(t)=Dist(X(:,t),Xstation)+sqrtm(R)*randn;%�۲�ֵ
end	

%EKF
Xekf=zeros(4,N);	
Xekf(:,1)=X(:,1);%Kalman�˲���״̬��ʼ��	
P0=eye(4);%���Э�������ĳ�ʼ��

for i=2:N
    Xn=F*Xekf(:,i-1);%һ��Ԥ��
    P1=F*P0*F'+G*Q*G';%Ԥ�����Э����
    dd=Dist(Xn,Xstation);%�۲�Ԥ��	
    %����ſ˱Ⱦ���H	
    H=[(Xn(1,1)-x0)/dd,0,(Xn(3,1)-y0)/dd,0];%̩��չ����һ�׽���
    K=P1*H'*inv(H*P1*H'+R);%��������������
    Xekf(:,i)=Xn+K*(Z(:,i)-dd);%״̬����
    P0=(eye(4)-K*H)*P1;%�˲����Э�������
end
	
%������	
for i=1:N
    Err_KalmanFilter(i)=Dist(X(:,i),Xekf(:,i));%
end 	
%��ͼ	
figure	
hold on;box on;	
plot(X(1,:),X(3,:),'-k');%��ʵ�켣	
plot(Xekf(1,:),Xekf(3,:),'-r');%��չKalman�˲��켣	
legend(' real trajectory','EKF trajectory');	
xlabel('X-axis  X/m');	
ylabel('Y-axis Y/m');
	
figure 	
hold on ;box on;	
plot(Err_KalmanFilter,'-ks','MarkerFace','r')
xlabel('Time /s');
ylabel('Position estimation bias   /m');

%�Ӻ��� ��ŷ�Ͼ���	
function d=Dist(X1,X2);	
    if length(X2)<=2 	
        d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);	
    else 	
        d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(3))^2);
    end
end
	
% %�Ӻ��� ��ŷ�Ͼ���
% function d=Dist(X1,X2)
%     if length(X2)<=2
%         d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(2))^2);
%     else	
%         d=sqrt((X1(1)-X2(1))^2+(X1(3)-X2(3))^2); 	
%     end
% end