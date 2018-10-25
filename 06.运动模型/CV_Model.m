%ת���������������˲�Ч��������Ӱ�죡��������

clc;
clear all;
close all;

T = 0.2;
SMP = 50;
F=[ 1 0 T 0; 
    0 1 0 T; 
    0 0 1 0; 
    0 0 0 1]; %״̬ת�ƾ���

dr1 = 0.4;
Q=[ dr1^2 0 0 0; 
    0 dr1^2 0 0; 
    0 0 dr1^2 0; 
    0 0 0 dr1^2]; %ת�������������

H=[ 1 0 0 0; 
    0 1 0 0]; %�������

dr2 = 2;
R=[ dr2^2 0 ; 
    0 dr2^2]; %���������������

P=[ 1 0 0 0; 
    0 1 0 0; 
    0 0 1 0; 
    0 0 0 1]; %Э�������

X(:,1)=[200;200;10;10]; %״̬
Z(:,1)= H*X(:,1)+dr2*randn(2,1);
for i=2:SMP  
    X(:,i)= F*X(:,i-1);                          %ģ�⣬����Ŀ���˶���ʵ�켣
    Z(:,i)= H*X(:,i)+dr2*randn(2,1);     %����ֵ��������
end

x = X(:,1);
for i=1:SMP
    [x,P] = kf(F,x,P,H,Z(:,i),Q,R);
    x_filter(:,i) = x;
end

 delta_r= sqrt( (x_filter(1,:)-X(1,:)).^2 + (x_filter(2,:)-X(2,:)).^2 );
 
figure;
plot(X(1,:),X(2,:),'-b.');hold on;
plot(x_filter(1,:),x_filter(2,:),'-r.');hold on;
plot(Z(1,:),Z(2,:),'*');hold on;grid on
xlabel('x(m)');ylabel('y(m)');
legend('ʵ��ֵ','�˲�ֵ','����ֵ');

figure;
plot(delta_r,'-b.');hold on;grid on
xlabel('x(m)');
ylabel('y(m)')


