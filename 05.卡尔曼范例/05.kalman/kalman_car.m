clc;
clear all;
close all;

X=[200;2000;15;0]; %״̬

F=[1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]; %״̬ת�ƾ���
% Q=[0.01 0 0 0; 0 0.0001 0 0; 0 0 0.0001 0; 0 0 0 0.0001]; %״̬ת��Э�������
Q = [0.01;0.01;0.0001;0.0001];
H=[1 0 0 0 ; 0 1 0 0; 0 0 1 0; 0 0 0 1]; %�۲����
% R=[1 0 0 0 ; 0 1 0 0; 0 0 0.0001 0; 0 0 0 0.0001]; %�۲���������
R = [1;1;0.0001;0.0001];
P=[1 0 0 0; 0 1 0 0 ; 0 0 1 0; 0 0 0 1]; %״̬Э�������

aa = zeros(4,100);
Z = zeros(4,100);
for i=2:100  
    X(:,i)= F*X(:,i-1)+Q; %ģ�⣬����Ŀ���˶���ʵ�켣
end

X1 = X(:,1);
aa(:,1) = X1;
for i=2:100
    Z(:,i)= H*X(:,i)+ R;     %����ֵ��������
    X_ = F*X1;
    P_ = F*P*F'+Q;
    K = P_*H'/(H*P_*H'+R);
    X1 = X_+K*(Z(:,i)-H*X_);
    P = (eye(4)-K*H)*P_;
  
    aa(:,i) = X1;
end

figure;
hold on;
plot(X(1,:),X(2,:)); 
hold on;
plot(aa(1,:),aa(2,:)); 

xlabel('x(m)');
ylabel('y(m)')
