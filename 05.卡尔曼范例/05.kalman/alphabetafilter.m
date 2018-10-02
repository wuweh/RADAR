%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ����������ϵ�������£���֤�������-���˲�����ϵ����۲��ź������йأ�����ο����ף�
%  ���룺��ʵ���ݣ�λ��+�ٶȣ�realval
%            �۲����ݣ�λ�ã�obseval �ӹ۲�����
%  ������������ݣ�λ��+�ٶȣ�estimval
%            У�����ֵ��λ��+�ٶȣ�reviseval
%  �ο����ף�������Ų��壬����Ȫ.������״����ݴ���̳�.���������ӹ�ҵ������.pp.60-66
%  ���Ĺ�/������Ϣ����ѧԺ/�������պ����ѧ
%  5/8/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;
clear all;

alpha=0.27;
beta=0.04229;

T=0.1;                      %�������/s
phi=[1,T;0,1];              %ת�ƾ���
K=[alpha;beta/T];           %�����������
H=[1,1];                    %��������

%% �������(����۲�ֵobseval)
t=0:T:100;                  %����ʱ��/s
LEN=length(t);              %�������ݳ���

xi=0;                       %λ�ù۲�������ֵ
zeta=40e3;                  %λ�ù۲���������
realval=zeros(2,LEN);                                 

realval(1,:)=0.3e3*t+100e3;         %λ����ʵֵ/m
realval(2,:)=0.3e3*ones(1,LEN);     %�ٶ���ʵֵ/m/s

origsite=100e3;                     %��ʼλ��/m
origvelo=0.3e3;                     %��ʼ�ٶ�/m/s

obseval=H*realval+(sqrt(zeta)*randn(1,LEN)+xi*ones(1,LEN)); %λ�ù۲�ֵ +����Ϊzeta��ֵΪ0�ĸ�˹����

figure(1);
subplot(221);plot(t,realval(1,:)/1e3);title('Ŀ����ʵ���');
xlabel('ʱ��/s');ylabel('����/km');grid on;

%% ���Ԥ����У��
estimval=zeros(2,LEN);       %���Ԥ��ֵ(����ֵ)
reviseval=zeros(2,LEN);      %���У��ֵ(���ֵ)

estimval(1,1)=origsite;        %��ʼ������ֵ��λ��/m��
estimval(2,1)=origvelo;       %��ʼ������ֵ���ٶ�/m/s��

reviseval(1,1)=origsite;        %��ʼ��У��ֵ��λ��/m��
reviseval(2,1)=origvelo;       %��ʼ��У��ֵ���ٶ�/m/s��

%%����������������һ���Ǿ��룻һ�����ٶ�
for n=1:1:LEN-1
    %��kʱ��У�����ֵԤ��k+1ʱ��״̬
    %estimval�ǿ�����ģ���е�Ԥ��ֵ
    %phi��ת�ƾ���
    estimval(:,n+1)=phi*reviseval(:,n);
    
    %��k+1ʱ�̻��ʵ�ʹ۲�ֵobseval���ڴ˻�����У�����ֵ
    %reviseval�ǿ�����ģ���е��������ֵ������KΪ����������ϵ�����ڦ�-���˲�ģ���У�Ϊһ����ֵ
    %H���������
    reviseval(:,n+1)=estimval(:,n+1)+K*(obseval(:,n+1)-H*estimval(:,n+1));
end

%% ������
subplot(222);
plot(t,estimval(1,:)/1e3);title('���ƹ��');
xlabel('ʱ��/s');ylabel('����/km');grid on;

subplot(223);
plot(t,reviseval(1,1:LEN)/1e3);title('У��(���)���');
xlabel('ʱ��/s');ylabel('����/km');grid on;

subplot(224);
plot(t,abs(estimval(1,1:LEN)-reviseval(1,1:LEN))/1e3);title('���');
xlabel('ʱ��/s');ylabel('���/km');grid on;

%λ��
figure;
plot(t,estimval(1,1:LEN)/1e3);hold on;
plot(t,reviseval(1,1:LEN)/1e3,'r');hold on;
plot(t,realval(1,1:LEN)/1e3,'k');hold off;legend('����ֵ','У��ֵ','�۲�ֵ');grid on;

%�ٶ�
figure;
plot(t,estimval(2,1:LEN));hold on;
plot(t,reviseval(2,1:LEN),'r');hold on;
plot(t,realval(2,1:LEN),'k');hold off;legend('����ֵ','У��ֵ','�۲�ֵ');