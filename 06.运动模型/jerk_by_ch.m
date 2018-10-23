
close all
clear all;
clc;
%������������������������������������������������������������������������������
T= 0.05;    %��������
SAMP=1000;   %��������
Pi=3.1415926536;
alpha=0.7;  %����Ƶ��
MK_num=50;
%������������������������������������������������������������������������������

X=zeros(8,SAMP);
X(:,1)=[10,2,1,0,30,20,1,0]';

W1= -20*Pi/180;
W2= 20*Pi/180; 

F1=[1,T,T^2/2,0,0,0,0,0;
    0,1,T,0,0,0,0,0;
    0,0,1,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,1,T,T^2/2,0;
    0,0,0,0,0,1,T,0;
    0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0];%����ֱ��ģ��

F2=[1,sin(W1*T)/W1,(1-cos(W1*T))/W1^2,0,0,0,0,0;
    0,cos(W1*T),sin(W1*T)/W1,0,0,0,0,0;
    0,-W1*sin(W1*T),cos(W1*T),0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,1,sin(W1*T)/W1,(1-cos(W1*T))/W1^2,0;
    0,0,0,0,0,cos(W1*T),sin(W1*T)/W1,0;
    0,0,0,0,0,-W1*sin(W1*T),cos(W1*T),0;
    0,0,0,0,0,0,0,0];

F3=[1,sin(W2*T)/W2,(1-cos(W2*T))/W2^2,0,0,0,0,0;
    0,cos(W2*T),sin(W2*T)/W2,0,0,0,0,0;
    0,-W2*sin(W2*T),cos(W2*T),0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,1,sin(W2*T)/W2,(1-cos(W2*T))/W2^2,0;
    0,0,0,0,0,cos(W2*T),sin(W2*T)/W2,0;
    0,0,0,0,0,-W2*sin(W2*T),cos(W2*T),0;
    0,0,0,0,0,0,0,0];%ת��ģ��

F4=[1,T,T^2/2,T^3/6,0,0,0,0;
    0,1,T,T^2/2,0,0,0,0;
    0,0,1,T,0,0,0,0;
    0,0,0,1,0,0,0,0
    0,0,0,0,1,T,T^2/2,T^3/6;
    0,0,0,0,0,1,T,T^2/2;
    0,0,0,0,0,0,1,T;
    0,0,0,0,0,0,0,1];%����ֱ��ģ��

for t=1:SAMP  
    if t==1 
        X(:,1)=[10,2,1,0,30,20,1,0]';
    elseif t>1&&t<100
        X(:,t)=F1*X(:,t-1);
    elseif t==100
        X(:,t)=F2*[X(1,t-1),X(2,t-1),-W1*X(6,t-1),0,X(5,t-1),X(6,t-1),W1*X(2,t-1),0]';
    elseif(t>100&&t<300)
        X(:,t)=F2*X(:,t-1);
    elseif (t==300)
        X(:,t)=F1*[X(1,t-1),X(2,t-1),0,0,X(5,t-1),X(6,t-1),0,0]';
    elseif (t>300&&t<400)
        X(:,t)=F1*X(:,t-1);
    elseif(t==400)
        X(:,t)=F3*[X(1,t-1),X(2,t-1),-W2*X(6,t-1),0,X(5,t-1),X(6,t-1),W2*X(2,t-1),0]';
    elseif(t>400&&t<700)
        X(:,t)=F3*X(:,t-1);
    elseif (t==700)
        X(:,t)=F1*[X(1,t-1),X(2,t-1),0,0,X(5,t-1),X(6,t-1),0,0]';
    elseif (t>700&&t<900)
        X(:,t)=F1*X(:,t-1);
    elseif(t==900)
        X(:,t)=F3*[X(1,t-1),X(2,t-1),-W2*X(6,t-1),0,X(5,t-1),X(6,t-1),W2*X(2,t-1),0]';
    else
        X(:,t)=F3*X(:,t-1);
    end
end

% t=1:SAMP;
% figure
% plot(X(1,t),X(5,t),'-b');
% axis equal;

%ת�ƾ���----------------------------------------------------------------------
PT=(2-2*alpha*T+alpha^2*T^2-2*exp(-alpha*T))/(2*alpha^3);
QT=(exp(-alpha*T)-1+alpha*T)/alpha^2;
RT=(1-exp(-alpha*T))/alpha;
ST=exp(-alpha*T);

FJT=[1,T,T^2/2,PT,0,0,0,0;
    0,1,T,QT,0,0,0,0;
    0,0,1,RT,0,0,0,0;
    0,0,0,ST,0,0,0,0;
    0,0,0,0,1,T,T^2/2,PT;
    0,0,0,0,0,1,T,QT;
    0,0,0,0,0,0,1,RT;
    0,0,0,0,0,0,0,ST];%F(T)Ϊģ�͵�ϵͳת�ƾ���
 
%״̬����Э�������------------------------------------------------------------
Q11=(alpha^5*T^5/10-alpha^4*T^4/2+4*alpha^3*T^3/3-2*alpha^2*T^2+2*alpha*T-3 ...
    +4*exp(-alpha*T)+2*alpha^2*T^2*exp(-alpha*T)-exp(-2*alpha*T))/(2*alpha^7);
Q12=(1-2*alpha*T+2*alpha^2*T^2-alpha^3*T^3+alpha^4*T^4/4+exp(-2*alpha*T)+ ...
    2*alpha*T*exp(-alpha*T)-2*exp(-alpha*T)-alpha^2*T^2*exp(-alpha*T))/(2*alpha^6);%��Q21��ͬ,����ͬ��
Q13=(2*alpha*T-alpha^2*T^2+alpha^3*T^3/3-3-exp(-2*-alpha*T)+4*exp(-alpha*T) ...
    +alpha^2*T^2*exp(-alpha*T))/(2*alpha^5);
Q14=(1+exp(-2*alpha*T)-2*exp(-alpha*T)-alpha^2*T^2*exp(-alpha*T))/(2*alpha^4);
Q21=Q12;
Q22=(1-exp(-2*alpha*T)+2*alpha^3*T^3/3+2*alpha*T-2*alpha^2*T^2-4*alpha*T*exp(-alpha*T))/(2*alpha^5);
Q23=(1+alpha^2*T^2-2*alpha*T+2*alpha*T*exp(-alpha*T)+exp(-2*alpha*T)-2*exp(-alpha*T))/(2*alpha^4);
Q24=(1-2*alpha*T*exp(-2*alpha*T)-exp(-2*alpha*T))/(2*alpha^3);
Q31=Q13;
Q32=Q23;
Q33=(4*exp(-alpha*T)-exp(-2*alpha*T)+2*alpha*T-3)/(2*alpha^3);
Q34=(1-2*exp(-alpha*T)+exp(-2*alpha*T))/(2*alpha^2);
Q41=Q14;
Q42=Q24;
Q43=Q34;
Q44=(1-exp(-2*alpha*T))/(2*alpha);
 
QJ=[Q11,Q12,Q13,Q14,0,0,0,0;
    Q21,Q22,Q23,Q24,0,0,0,0;
    Q31,Q32,Q33,Q34,0,0,0,0;
    Q41,Q42,Q43,Q44,0,0,0,0;
    0,0,0,0,Q11,Q12,Q13,Q14;
    0,0,0,0,Q21,Q22,Q23,Q24;
    0,0,0,0,Q31,Q32,Q33,Q34;
    0,0,0,0,Q41,Q42,Q43,Q44;];

HJ=[1,0,0,0,0,0,0,0;
    0,0,0,0,1,0,0,0];

RX=2; % X����۲�����׼��
RY=2; % Y����۲�����׼��
R=[RX^2,0;
   0,RY^2];%�۲�����Э�������


%����������Ĺ켣=�۲�ֵ ������������������������������������������������������
noise_x=RX*randn(1,SAMP);%randn�����ľ�ֵΪ0,1Ϊ��׼�SAMP�������ֵ
Z_X=X(1,:)+noise_x; %X���������
noise_y=RY*randn(1,SAMP);
Z_Y=X(5,:)+noise_y; %Y���������
Z=[Z_X;Z_Y];

XJk=zeros(8,SAMP);
PJk=zeros(8,8,SAMP);

XJk_Initial=[X(1,1),X(2,1),X(3,1),0,X(5,1),X(6,1),X(7,1),0]';

PJ_Initial=[   RX^2,0,0,0,0,0,0,0;
               0,10,0,0,0,0,0,0;
               0,0,10,0,0,0,0,0;
               0,0,0,1,0,0,0,0;
               0,0,0,0,RY^2,0,0,0;
               0,0,0,0,0,10,0,0;
               0,0,0,0,0,0,10,0;
               0,0,0,0,0,0,0,1];%�������Э�������
           
for i=1:MK_num
        PJk(:,:,1)=PJ_Initial;
        XJk(:,1)=XJk_Initial;
        
    for t=1:1:SAMP           
        [XJk(:,t),PJk(:,:,t)] = kf(FJT,XJk(:,t),PJk(:,:,t),HJ,Z(:,t),QJ,R);     %����jerkģ��            
%        [XJk(:,t),PJk(:,:,t)] = kf(F,a,PJk(:,:,t),HJ,Z(:,t),QJ_2,R);     %�����ȼ���ģ�� 
        a = XJk(:,t);
        b = PJk(:,:,t);
    end
     XJK_all(:,:,i)=XJk;
end

XJk=sum(XJK_all,3)/MK_num;

% % D=sqrt((Xk(1,:)-X(1,:)).^2+(Xk(6,:)-X(5,:)).^2);%���������
 DJ=sqrt((XJk(1,:)-X(1,:)).^2+(XJk(5,:)-X(5,:)).^2);%���������
% 
% % V=sqrt((Xk(2,:)-X(2,:)).^2+(Xk(7,:)-X(6,:)).^2);%���������
 VJ=sqrt((XJk(2,:)-X(2,:)).^2+(XJk(6,:)-X(6,:)).^2);%���������
% 
% % A=sqrt((Xk(3,:)-X(3,:)).^2+(Xk(8,:)-X(7,:)).^2);%���������
 AJ=sqrt((XJk(3,:)-X(3,:)).^2+(XJk(7,:)-X(7,:)).^2);%���������
% 
% % a=sqrt((Xk(4,:)-X(4,:)).^2+(Xk(9,:)-X(8,:)).^2);%���������
 aJ=sqrt((XJk(4,:)-X(4,:)).^2+(XJk(8,:)-X(8,:)).^2);%���������

t=1:SAMP;

figure
plot(Z_X,Z_Y,'r.');
hold on;
plot(X(1,t),X(5,t),'-g.');
hold on;
plot(XJk(1,t),XJk(5,t),'-b.');
xlabel('x(t)(m)'),ylabel('y(t)(m)'); 
legend('����λ��','ʵ��λ��','����λ��');
title('Jerk�㷨');
axis equal;%�����᳤�ȵ�λ������

% figure
% % plot(D(1:599),'b-');
% hold on;
% plot(DJ(1:599),'r-');
% xlabel('t'),ylabel('���������(m)'); 
% legend('��-Jerk�㷨','Jerk�㷨');
% title('λ�ù������');
% 
% figure
% % plot(V(1:599),'b-');
% hold on;
% plot(VJ(1:599),'r-');
% xlabel('t'),ylabel('V(m/s)'); 
% legend('��-Jerk�㷨','Jerk�㷨');
% title('�ٶȹ������');
% 
% figure
% % plot(A(1:599),'b-');
% hold on;
% plot(AJ(1:599),'r-');
% xlabel('t'),ylabel('V(m^2/s)'); 
% legend('��-Jerk�㷨','Jerk�㷨');
% title('���ٶȹ������');
% 
% figure
% % plot(a(1:599),'b-');
% hold on;
% plot(aJ(1:599),'r-');
% xlabel('t'),ylabel('V(m^2/s)'); 
% legend('��-Jerk�㷨','Jerk�㷨');
% title('�����������');
