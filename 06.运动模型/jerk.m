close all
clear all;
clc;
%———————————————————————————————————————
T= 0.5;%采样周期
SAMP=600;%采样点数
Pi=3.1415926536;
alpha=0.5;%急动频率
MK_num=200;
%———————————————————————————————————————

%原始轨迹,前35秒作大约为1.88马赫直线运动，35-94秒过载1.22g，94-106秒1.88马赫直线运动，106-136秒过载3.7g,136-178秒1.88马赫直线运动，178-190秒过载
%6g,190-200秒1.88马赫直线运动.
X=zeros(8,SAMP);
X(:,1)=[30000,300,0,0,60000,400,0,0]';%给X状态赋初值；：作用，矩阵X中，取所有行第一列的元素

W1= 1.4*3.1415926536/180;
W2= -4.2*3.1415926536/180;%各转角的角速度

F1=[1,T,T^2/2,0,0,0,0,0;
    0,1,T,0,0,0,0,0;
    0,0,1,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,1,T,T^2/2,0;
    0,0,0,0,0,1,T,0;
    0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0];%匀速直线模型

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
    0,0,0,0,0,0,0,0];%转弯模型

F4=[1,T,T^2/2,T^3/6,0,0,0,0;
    0,1,T,T^2/2,0,0,0,0;
    0,0,1,T,0,0,0,0;
    0,0,0,1,0,0,0,0
    0,0,0,0,1,T,T^2/2,T^3/6;
    0,0,0,0,0,1,T,T^2/2;
    0,0,0,0,0,0,1,T;
    0,0,0,0,0,0,0,1];%匀速直线模型

for t=1:SAMP  %从时刻2开始，到采样结束为止
    if t==1 
        X(:,1)=[30000,300,0,0,60000,400,0,0]';
    elseif t>1&&t<150
        X(:,t)=F1*X(:,t-1);
    elseif t==150
        X(:,t)=F2*[X(1,t-1),X(2,t-1),-W1*X(6,t-1),0,X(5,t-1),X(6,t-1),W1*X(2,t-1),0]';
    elseif(t>150&&t<300)
        X(:,t)=F2*X(:,t-1);
    elseif (t==300)
        X(:,t)=F1*[X(1,t-1),X(2,t-1),0,0,X(5,t-1),X(6,t-1),0,0]';
    elseif (t>300&&t<450)
        X(:,t)=F1*X(:,t-1);
    elseif(t==450)
        X(:,t)=F2*[X(1,t-1),X(2,t-1),-W2*X(6,t-1),0,X(5,t-1),X(6,t-1),W2*X(2,t-1),0]';
    elseif(t>450&&t<540)
        X(:,t)=F4*X(:,t-1);
    else
        X(:,t)=F1*X(:,t-1);
    end
end

%改进的当前模型——————分别对X和Y方向进行建模—————————————————————
XT=zeros(10,SAMP);%系统模型
X(:,1)=[30000,300,0,0,60000,400,0,0]';

XT(1,SAMP)=X(1,1);%给模型赋X方向的位移初值
XT(6,SAMP)=X(5,1);
XT(2,SAMP)=X(2,1);%给模型赋X方向的速度初值
XT(7,SAMP)=X(6,1);
XT(3,SAMP)=X(3,1);%给模型赋X方向的加速度初值
XT(8,SAMP)=X(7,1);
% XT(4,SAMP)=0;
% XT(9,SAMP)=0;
XT(5,SAMP)=alpha;
XT(10,SAMP)=alpha;
%转移矩阵----------------------------------------------------------------------
PT=(2-2*alpha*T+alpha^2*T^2-2*exp(-alpha*T))/(2*alpha^3);
QT=(exp(-alpha*T)-1+alpha*T)/alpha^2;
RT=(1-exp(-alpha*T))/alpha;
ST=exp(-alpha*T);

U1=(2*T-alpha*T^2+alpha^3*T^3/3-2*(1-exp(-alpha*T))/alpha)/(2*alpha^2);
U2=(-T+alpha*T^2/2+(1-exp(-alpha*T))/alpha)/alpha;
U3=T-(1-exp(-alpha*T))/alpha;
U4=1-exp(-alpha*T);
U=[U1,U2,U3,U4,0,U1,U2,U3,U4,0]';

F=[1,T,T^2/2,T^3/6,0,0,0,0,0,0;
   0,1,T,T^2/2,0,0,0,0,0,0;
   0,0,1,T,0,0,0,0,0,0;
   0,0,0,1,0,0,0,0,0,0;
   0,0,0,0,1,0,0,0,0,0;
   0,0,0,0,0,1,T,T^2/2,T^3/6,0;
   0,0,0,0,0,0,1,T,T^2/2,0;
   0,0,0,0,0,0,0,1,T,0;
   0,0,0,0,0,0,0,0,1,0;
   0,0,0,0,0,0,0,0,0,1];
FT=[1,T,T^2/2,PT,0,0,0,0,0,0;
    0,1,T,QT,0,0,0,0,0,0;
    0,0,1,RT,0,0,0,0,0,0;
    0,0,0,ST,0,0,0,0,0,0;
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,T,T^2/2,PT,0;
    0,0,0,0,0,0,1,T,QT,0;
    0,0,0,0,0,0,0,1,RT,0;
    0,0,0,0,0,0,0,0,ST,0;
    0,0,0,0,0,0,0,0,0,1];%F(T)为模型的系统转移矩阵

FJT=[1,T,T^2/2,PT,0,0,0,0;
    0,1,T,QT,0,0,0,0;
    0,0,1,RT,0,0,0,0;
    0,0,0,ST,0,0,0,0;
    0,0,0,0,1,T,T^2/2,PT;
    0,0,0,0,0,1,T,QT;
    0,0,0,0,0,0,1,RT;
    0,0,0,0,0,0,0,ST];%F(T)为模型的系统转移矩阵

%状态噪声协方差矩阵------------------------------------------------------------
Q11=(alpha^5*T^5/10-alpha^4*T^4/2+4*alpha^3*T^3/3-2*alpha^2*T^2+2*alpha*T-3+4*exp(-alpha*T)+2*alpha^2*T^2*exp(-alpha*T)-exp(-2*alpha*T))/(2*alpha^7);
Q12=(1-2*alpha*T+2*alpha^2*T^2-alpha^3*T^3+alpha^4*T^4/4+exp(-2*alpha*T)+2*alpha*T*exp(-alpha*T)-2*exp(-alpha*T)-alpha^2*T^2*exp(-alpha*T))/(2*alpha^6);%与Q21相同,以下同理
Q13=(2*alpha*T-alpha^2*T^2+alpha^3*T^3/3-3-exp(-2*-alpha*T)+4*exp(-alpha*T)+alpha^2*T^2*exp(-alpha*T))/(2*alpha^5);
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

%Sigmaalpha的取值----------------------------------------------------
Sigmaj=0.2;
SigmaW=sqrt(2*alpha*Sigmaj^2);
SigmaEpsilong=0.7;
Q=[Q11*SigmaW^2,Q12*SigmaW^2,Q13*SigmaW^2,Q14*SigmaW^2,0,0,0,0,0,0;
     Q21*SigmaW^2,Q22*SigmaW^2,Q23*SigmaW^2,Q24*SigmaW^2,0,0,0,0,0,0; 
     Q31*SigmaW^2,Q32*SigmaW^2,Q33*SigmaW^2,Q34*SigmaW^2,0,0,0,0,0,0;
     Q41*SigmaW^2,Q42*SigmaW^2,Q43*SigmaW^2,Q44*SigmaW^2,0,0,0,0,0,0;
     0,0,0,0,T*SigmaEpsilong^2,0,0,0,0,0;
     0,0,0,0,0,Q11*SigmaW^2,Q12*SigmaW^2,Q13*SigmaW^2,Q14*SigmaW^2,0;
     0,0,0,0,0,Q21*SigmaW^2,Q22*SigmaW^2,Q23*SigmaW^2,Q24*SigmaW^2,0;
     0,0,0,0,0,Q31*SigmaW^2,Q32*SigmaW^2,Q33*SigmaW^2,Q34*SigmaW^2,0;
     0,0,0,0,0,Q41*SigmaW^2,Q42*SigmaW^2,Q43*SigmaW^2,Q44*SigmaW^2,0;
     0,0,0,0,0,0,0,0,0,T*SigmaEpsilong^2];%状态噪声协方差阵
QJ=[Q11,Q12,Q13,Q14,0,0,0,0;
    Q21,Q22,Q23,Q24,0,0,0,0;
    Q31,Q32,Q33,Q34,0,0,0,0;
    Q41,Q42,Q43,Q44,0,0,0,0;
    0,0,0,0,Q11,Q12,Q13,Q14;
    0,0,0,0,Q21,Q22,Q23,Q24;
    0,0,0,0,Q31,Q32,Q33,Q34;
    0,0,0,0,Q41,Q42,Q43,Q44;];

H=[1,0,0,0,0,0,0,0,0,0;
   0,0,0,0,0,1,0,0,0,0];%观测矩阵
HJ=[1,0,0,0,0,0,0,0;
    0,0,0,0,1,0,0,0];

RX=100; % X方向观测误差标准差
RY=100; % Y方向观测误差标准差
R=[RX^2,0;
   0,RY^2];%观测噪声协方差矩阵


%叠加噪声后的轨迹=观测值 ———————————————————————————
noise_x=RX*randn(1,SAMP);%randn函数的均值为0,1为标准差，SAMP处代表均值
Z_X=X(1,:)+noise_x; %X方向加噪声
noise_y=RY*randn(1,SAMP);
Z_Y=X(5,:)+noise_y; %Y方向加噪声
Z=[Z_X;Z_Y];
ZZ=Z(:,t);

XK_all=zeros(10,SAMP,MK_num);
XJK_all=zeros(8,SAMP,MK_num);
Xk=zeros(10,SAMP);
XJk=zeros(8,SAMP);
Pk=zeros(10,10,SAMP);
PJk=zeros(8,8,SAMP);
K=zeros(10,2,SAMP);%存放增益矩阵

Xk_Initial=[X(1,1),X(2,1),X(3,1),0,alpha,X(5,1),X(6,1),X(7,1),0,alpha]';
XJk_Initial=[X(1,1),X(2,1),X(3,1),0,X(5,1),X(6,1),X(7,1),0]';

PJ_Initial=[RX^2,0,0,0,0,0,0,0;
   0,1000,0,0,0,0,0,0;
   0,0,100,0,0,0,0,0;
   0,0,0,10,0,0,0,0;
   0,0,0,0,RY^2,0,0,0;
   0,0,0,0,0,1000,0,0;
   0,0,0,0,0,0,100,0;
   0,0,0,0,0,0,0,10];%估计误差协方差矩阵

P_Initial=[RX^2,0,0,0,0,0,0,0,0,0;
   0,1000,0,0,0,0,0,0,0,0;
   0,0,100,0,0,0,0,0,0,0;
   0,0,0,10,0,0,0,0,0,0;
   0,0,0,0,10,0,0,0,0,0;
   0,0,0,0,0,RY^2,0,0,0,0;
   0,0,0,0,0,0,1000,0,0,0;
   0,0,0,0,0,0,0,100,0,0;
   0,0,0,0,0,0,0,0,10,0;
   0,0,0,0,0,0,0,0,0,10];%估计误差协方差矩阵

for i=1:MK_num

        Pk(:,:,1)=P_Initial;
        Xk(:,1)=Xk_Initial;
        PJk(:,:,1)=PJ_Initial;
        XJk(:,1)=XJk_Initial;
        P0=P_Initial;
        PJ0=PJ_Initial;
        Xk0=Xk_Initial;
        XkJ0=XJk_Initial;
for t=2:1:SAMP-1
    

        
     Xk(:,t)=F*Xk0;%计算一步预测估计值
     Pk(:,:,t)=FT*P0*FT'+Q;%计算预测滤波协方差矩阵
     K(:,:,t)=Pk(:,:,t)*H'*inv(H*Pk(:,:,t)*H'+R);%计算增益矩阵
     Xk(:,t+1)=Xk(:,t)+K(:,:,t)*(Z(:,t)-H*Xk(:,t));%计算估计值
     Pk(:,:,t+1)=(eye(10)-K(:,:,t)*H)*Pk(:,:,t);%计算协方差更新方程
     P0=Pk(:,:,t+1);%估计误差协方差阵更新
     Xk0=Xk(:,t+1);%估计矩阵更新
     
     XJk(:,t)=FJT*XkJ0;%计算一步预测估计值
     PJk(:,:,t)=FJT*PJ0*FJT'+QJ;%计算预测滤波协方差矩阵
     KJ(:,:,t)=PJk(:,:,t)*HJ'*inv(HJ*PJk(:,:,t)*HJ'+R);%计算增益矩阵
     XJk(:,t+1)=XJk(:,t)+KJ(:,:,t)*(Z(:,t)-HJ*XJk(:,t));%计算估计值
     PJk(:,:,t+1)=(eye(8)-KJ(:,:,t)*HJ)*PJk(:,:,t);%计算协方差更新方程
     PJ0=PJk(:,:,t+1);%估计误差协方差阵更新
     XkJ0=XJk(:,t+1);%估计矩阵更
     
     

end
    XK_all(:,:,i)=Xk;
     XJK_all(:,:,i)=XJk;
end
Xk=sum(XK_all,3)/MK_num;
XJk=sum(XJK_all,3)/MK_num;

D=sqrt((Xk(1,:)-X(1,:)).^2+(Xk(6,:)-X(5,:)).^2);%均方根误差
DJ=sqrt((XJk(1,:)-X(1,:)).^2+(XJk(5,:)-X(5,:)).^2);%均方根误差

V=sqrt((Xk(2,:)-X(2,:)).^2+(Xk(7,:)-X(6,:)).^2);%均方根误差
VJ=sqrt((XJk(2,:)-X(2,:)).^2+(XJk(6,:)-X(6,:)).^2);%均方根误差

A=sqrt((Xk(3,:)-X(3,:)).^2+(Xk(8,:)-X(7,:)).^2);%均方根误差
AJ=sqrt((XJk(3,:)-X(3,:)).^2+(XJk(7,:)-X(7,:)).^2);%均方根误差

a=sqrt((Xk(4,:)-X(4,:)).^2+(Xk(9,:)-X(8,:)).^2);%均方根误差
aJ=sqrt((XJk(4,:)-X(4,:)).^2+(XJk(8,:)-X(8,:)).^2);%均方根误差

t=1:SAMP;
figure
plot(Z_X,Z_Y,'-r');
hold on;
plot(X(1,t),X(5,t),'-g');
hold on;
plot(Xk(1,t),Xk(6,t),'-b');
xlabel('x(t)(m)'),ylabel('y(t)(m)'); 
legend('测量位置','实际位置','估计位置');
title('α-Jerk算法');
axis equal;%坐标轴长度单位设成相等

figure
plot(Z_X,Z_Y,'-r');
hold on;
plot(X(1,t),X(5,t),'-g');
hold on;
plot(XJk(1,t),XJk(5,t),'-b');
xlabel('x(t)(m)'),ylabel('y(t)(m)'); 
legend('测量位置','实际位置','估计位置');
title('Jerk算法');
axis equal;%坐标轴长度单位设成相等

figure
plot(D(1:599),'b-');
hold on;
plot(DJ(1:599),'r-');
xlabel('t'),ylabel('均方根误差(m)'); 
legend('α-Jerk算法','Jerk算法');
title('位置估计误差');

figure
plot(V(1:599),'b-');
hold on;
plot(VJ(1:599),'r-');
xlabel('t'),ylabel('V(m/s)'); 
legend('α-Jerk算法','Jerk算法');
title('速度估计误差');

figure
plot(A(1:599),'b-');
hold on;
plot(AJ(1:599),'r-');
xlabel('t'),ylabel('V(m^2/s)'); 
legend('α-Jerk算法','Jerk算法');
title('加速度估计误差');

figure
plot(a(1:599),'b-');
hold on;
plot(aJ(1:599),'r-');
xlabel('t'),ylabel('V(m^2/s)'); 
legend('α-Jerk算法','Jerk算法');
title('急动估计误差');
