%
% ����IMM�㷨��Ŀ�����
clc;clear all;close all;

%===============================
%����ģ��
%===============================
% �������
simTime=200;      %�����������
T=1;                     %����ʱ��
w2=7*2*pi/360;     %ģ��2ת����3��
w3=-6*2*pi/360;    %ģ��3ת����-3��
H=[1,0,0,0;0,0,1,0];                      %ģ���������
G=[T^2/2,0;T,0;0,T^2/2;0,T];              %ģ�͹���������Ȩ����
r=3^2;                                 %20 2000
R=[r,0;0,r];                            %ģ����������Э�������
Q=[2^2,0;0,2^2];                                  %ģ�͹�������Э�������

F1=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];     %ģ��1״̬ת�ƾ���

F2=[1,sin(w2*T)/w2,0,(cos(w2*T)-1)/w2;
    0,cos(w2*T),0,sin(w2*T);
    0,(1-cos(w2*T))/w2,1,sin(w2*T)/w2;
    0,-sin(w2*T),0,cos(w2*T)];            %ģ��2״̬ת�ƾ��� ��ת��

F3=[1,sin(w3*T)/w3,0,(cos(w3*T)-1)/w3;
    0,cos(w3*T),0,sin(w3*T);
    0,(1-cos(w3*T))/w3,1,sin(w3*T)/w3;
    0,-sin(w3*T),0,cos(w3*T)];            %ģ��3״̬ת�ƾ��� ��ת��

x0=[100,20,100,20]';  % ��ʼ״̬

% ������������
x = zeros(4,simTime);
z = zeros(2,simTime);         %��������������
z_true = zeros(2,simTime);    %��ֵ����
measureNoise = sqrt(R)*randn(2,simTime);  %������������
x(:,1)=x0;
z(:,1)=H*x(:,1)+measureNoise(:,1);
z_true(:,1)=H*x(:,1);
for a=2:simTime
    if (a>=50)&&(a<=120) 
        x(:,a)=F2*x(:,a-1);      %20--40s��ת 
    elseif (a>=140)&&(a<=200) 
       x(:,a)=F3*x(:,a-1);        %60--80s��ת 
    else
        x(:,a)=F1*x(:,a-1);      %����ֱ���˶�
    end
z(:,a)=H*x(:,a)+measureNoise(:,a);
z_true(:,a)=H*x(:,a);
end
%===============================
%     IMM
%===============================

%��ʼ��
x1_IMM = zeros(4,1);      %ģ��1IMM�㷨״̬����ֵ
x2_IMM = zeros(4,1);      %ģ��2IMM�㷨״̬����ֵ
x3_IMM = zeros(4,1);      %ģ��3IMM�㷨״̬����ֵ
x_pro_IMM = zeros(4,simTime);   %IMM�㷨ģ���ۺ�״̬����ֵ
P_IMM=zeros(4,4,simTime);       %IMM�㷨ģ���ۺ�״̬Э�������
P1_IMM=zeros(4,4);
P2_IMM=zeros(4,4);
P3_IMM=zeros(4,4);              %IMM�㷨��ģ��Э�������
r1_IMM=zeros(2,1);
r2_IMM=zeros(2,1);
r3_IMM=zeros(2,1);
S1_IMM=zeros(2,2);
S2_IMM=zeros(2,2);
S3_IMM=zeros(2,2);

%��ʼ��
x_pro_IMM(:,1)=x0;

pij=[0.95,0.025,0.025;
     0.025,0.95,0.025;
     0.025,0.025,0.95];    %ģ��ת�Ƹ��ʾ���
% pij=[0.7,0.2,0.1;
%      0.2,0.7,0.1;
%      0.1,0.2,0.7];    %ģ��ת�Ƹ��ʾ���

u_IMM=zeros(3,simTime);
u_IMM(:,1)=[1/3,1/3,1/3]';  %IMM�㷨ģ�͸���
x1_IMM=x0;x2_IMM=x0;x3_IMM=x0;  %IMM�㷨��ģ�ͳ�ʼ״̬

P0=diag([10,10,10,10]);  %��ʼ״̬Э�������
P1_IMM=P0;P2_IMM=P0;P3_IMM=P0;
P_IMM(:,:,1)=P0;

x4 = x1_IMM;
x4_rec(:,1) = x4;
P4 = P0;
%main loop
for t=1:simTime-1
    %����tʱ�̴���ÿ��ģ�͵ĸ���
    c_j=pij'*u_IMM(:,t); 
    
    %����ÿ��ģ�͵���������
    ui1=(1/c_j(1))*pij(:,1).*u_IMM(:,t);
    ui2=(1/c_j(2))*pij(:,2).*u_IMM(:,t);
    ui3=(1/c_j(3))*pij(:,3).*u_IMM(:,t);    
    
    % �����ģ���˲���ʼ������
    x01=x1_IMM*ui1(1)+x2_IMM*ui1(2)+x3_IMM*ui1(3);
    x02=x1_IMM*ui2(1)+x2_IMM*ui2(2)+x3_IMM*ui2(3);
    x03=x1_IMM*ui3(1)+x2_IMM*ui3(2)+x3_IMM*ui3(3);   %��ģ���˲���ʼ״̬
    
    P01=(P1_IMM+(x1_IMM-x01)*(x1_IMM-x01)')*ui1(1)+...
        (P2_IMM+(x2_IMM-x01)*(x2_IMM-x01)')*ui1(2)+...
        (P3_IMM+(x3_IMM-x01)*(x3_IMM-x01)')*ui1(3);
    P02=(P1_IMM+(x1_IMM-x02)*(x1_IMM-x02)')*ui2(1)+...
        (P2_IMM+(x2_IMM-x02)*(x2_IMM-x02)')*ui2(2)+...
        (P3_IMM+(x3_IMM-x02)*(x3_IMM-x02)')*ui2(3);
    P03=(P1_IMM+(x1_IMM-x03)*(x1_IMM-x03)')*ui3(1)+...
        (P2_IMM+(x2_IMM-x03)*(x2_IMM-x03)')*ui3(2)+...
        (P3_IMM+(x3_IMM-x03)*(x3_IMM-x03)')*ui3(3);  %��ģ���˲���ʼ״̬Э�������
    
    %�ڶ���--�������˲�
    %ģ��1,2,3��׼�������˲�
     [x1_IMM,P1_IMM,r1_IMM,S1_IMM]=Kalman(x01,P01,z(:,t+1),F1,G,Q,H,R);
     [x2_IMM,P2_IMM,r2_IMM,S2_IMM]=Kalman(x02,P02,z(:,t+1),F2,G,Q,H,R);
     [x3_IMM,P3_IMM,r3_IMM,S3_IMM]=Kalman(x03,P03,z(:,t+1),F3,G,Q,H,R);
     
     %ģ��1,2,3-->>UKF�˲�
%      [x1_IMM,P1_IMM,r1_IMM,S1_IMM] = ukf(F1,x01,P01,H,z(:,t+1),G*Q*G',R);
%      [x2_IMM,P2_IMM,r2_IMM,S2_IMM] = ukf(F2,x02,P02,H,z(:,t+1),G*Q*G',R);
%      [x3_IMM,P3_IMM,r3_IMM,S3_IMM] = ukf(F3,x03,P03,H,z(:,t+1),G*Q*G',R);
    
    %������--ģ�͸��ʸ���
    [u_IMM(:,t+1)]=Model_P_up(r1_IMM,r2_IMM,r3_IMM,S1_IMM,S2_IMM,S3_IMM,c_j);
    
    %���Ĳ�--ģ���ۺ�
    [x_pro_IMM(:,t+1),P_IMM(:,:,t+1)]=Model_mix(x1_IMM,x2_IMM,x3_IMM,P1_IMM,P2_IMM,P3_IMM,u_IMM(:,t));
    
    %��IMM�㷨
    [x4,P4,r3_IMM,S3_IMM] = Kalman(x4,P4,z(:,t+1),F1,G,Q,H,R);
    x4_rec(:,t+1) = x4;
end


%===============================
%��ͼ
figure(1)
plot(z_true(1,:),z_true(2,:)); grid on; hold on
plot(x_pro_IMM(1,:),x_pro_IMM(3,:),'r');
plot(z(1,:),z(2,:),'*'); 
plot(x4_rec(1,:),x4_rec(3,:),'b'); 
hold off
title('Ŀ���˶��켣');
xlabel('x/m'); ylabel('y/m');
legend('��ʵֵ','�˲�ֵ','����ֵ');text(z(1,1)+500,z(2,1),'t=1');

% λ�����
figure(2)
subplot(2,1,1);
t=1:simTime;
plot(t,abs(x_pro_IMM(1,t)-x(1,t)),'LineWidth',1);hold on
plot(t,abs(x4_rec(1,t)-x(1,t)),'LineWidth',1);hold on;grid on;
title('x����λ�ø������');
xlabel('t/s'); ylabel('x-error/m');
legend('IMM','CT');

subplot(2,1,2);
t=1:simTime;
plot(t,abs(x_pro_IMM(3,t)-x(3,t)),'LineWidth',1);hold on;
plot(t,abs(x4_rec(3,t)-x(3,t)),'LineWidth',1);hold on;grid on;
title('y����λ�ø������');
xlabel('t/s'); ylabel('y-error/m');
legend('IMM','CT');

% % �ٶ����
% figure(3)
% subplot(2,1,1);
% t=1:simTime;
% plot(t,abs(x_pro_IMM(2,t)-x(2,t)),'LineWidth',2);grid on
% title('x�����ٶȸ������');
% xlabel('t/s'); ylabel('vx-error/m');
% 
% subplot(2,1,2);
% t=1:simTime;
% plot(t,abs(x_pro_IMM(4,t)-x(4,t)),'LineWidth',2);grid on
% title('y�����ٶȸ������');
% xlabel('t/s'); ylabel('vy-error/m');
% 
% ģ�͸���
t=1:simTime;
figure(4)
plot(t,u_IMM(1,t),'k.-',t,u_IMM(2,t),'r.-',t,u_IMM(3,t),'b.-');grid on
title('IMM�㷨ģ�͸�������');
axis([0 simTime-1 0 1])
xlabel('t/s'); ylabel('ģ�͸���');
legend('ģ��1','ģ��2','ģ��3');


%����ģ�͸��ʸ��º���
function [u]=Model_P_up(r1,r2,r3,S1,S2,S3,c_j)
%ģ�͸��ʸ��º���

%u  ģ�͸���    
%r1 ģ��1Ԥ�����
%r2 ģ��2Ԥ�����
%r3 ģ��3Ԥ�����
%S1 ģ��1Ԥ�����Э�������
%S2 ģ��2Ԥ�����Э�������
%S3 ģ��3Ԥ�����Э�������
%c_j  ģ�ͻ�ϸ���

%������Ȼ����
Lfun1=(1/sqrt(abs(2*pi*(det(S1)))))*exp((-1/2)*(r1'*inv(S1)*r1));   %Lfun1=1/(r1'*inv(S1)*r1);
Lfun2=(1/sqrt(abs(2*pi*(det(S2)))))*exp((-1/2)*(r2'*inv(S2)*r2));   %Lfun2=1/(r2'*inv(S2)*r2);
Lfun3=(1/sqrt(abs(2*pi*(det(S3)))))*exp((-1/2)*(r3'*inv(S3)*r3));    %Lfun3=1/(r3'*inv(S3)*r3);
% ����ģ�͸��¸��ʣ���������һʱ��ģ�͸���
c=[Lfun1,Lfun2,Lfun3]'.*c_j;
% �ٹ�һ��
u=[Lfun1,Lfun2,Lfun3]'.*c_j*(1/sum(c));

end

%����״̬��Ϻ���
function [x_pro,P]=Model_mix(x1,x2,x3,P1,P2,P3,u)
% ģ���ۺϺ���
%x_pro  ״̬�ۺ�ֵ
%P      �ۺ�Э�������
%x1     ģ��1״̬����ֵ
%x2     ģ��2״̬����ֵ
%x3     ģ��3״̬����ֵ
%P1     ģ��1״̬����Э�������
%P2     ģ��2״̬����Э�������
%P3     ģ��3״̬����Э�������
%u      ģ��ת������

% �����ʼ�ȨOK
x_pro=x1*u(1)+x2*u(2)+x3*u(3);
P=(P1+(x1-x_pro)*(x1-x_pro)')*u(1)+...
  (P2+(x2-x_pro)*(x2-x_pro)')*u(2)+...
  (P3+(x3-x_pro)*(x3-x_pro)')*u(3);
end


%�ġ�Kalman�˲�����
function [X,P,e,S]=Kalman(X_Forward,P_Forward,Z,A,G,Q,H,R)
%�������˲�2012.2.27   IMMר�ã��������в�ͬ
%����˵��
%       Z--�۲�����ʸ��
%       A--ϵͳģ��״̬����
%       G--ϵͳģ������ϵ������
%       Q--ϵͳģ����������
%       H--����ϵ������
%       R--����ģ������Э����
%       X_Forward--ǰ�ι���״̬ʸ��
%       P_Forward--ǰ�ι���״̬Э�������

%       X--�������״̬ʸ��
%       P--�������״̬Э�������
%       e--�в�
%       S--�в�Э�������

    % Ԥ��
    X_Pre=A*X_Forward;
    P_Pre=A*P_Forward*A'+G*Q*G';
    e = Z-H*(A*X_Forward); %�в�
    S =H*P_Pre*H'+R;  %�в�Э�������
    % �������
    K=P_Pre*H'*inv(S);

    % �����˲�ֵ�����Э������
    X=X_Pre+K*e;

    M=K*H;
    n=size(M);
    I=eye(n);
    P=(I-K*H)*P_Pre*(I-K*H)'+ K*R*K';
end


