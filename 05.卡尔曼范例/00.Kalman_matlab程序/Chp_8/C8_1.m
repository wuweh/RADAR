%C8_1.m
clc
clear
summ=0;
N=10;  %%%%计算10次，在求平均，需要在GUI上设置个可以设置数字的文板框，再把这个数据提出来应用到程序中
covvv=[];   %%%估计方差
for n=1:N
%%%%%不规则采样跟踪开始
%%%%%读取数据
load RFIDm5 
TT=[];TTT=[];  %%%为记录下面计算得出的量，需要先设置两个空变量
%%%%%%%%%%%%%计算得到不规则采样数据的周期
for i=1:length(dm(1,:))
if i==1
        T=ts(1);tt=0;
         TT=[TT T];TTT=[TTT tt];
    elseif i==2;
         T=ts(i)-ts(i-1);tt=0;
         TT=[TT T];TTT=[TTT tt];
    else
        T=ts(i)-ts(i-1);
         TT=[TT T];
         tt=TT(i)-TT(i-1);TTT=[TTT tt]; 
end
end
 
%%%%%%%%%%%选择模型参数
 ax=1/20;
 xamax=3;
 qqx=(xamax)^2*(4-pi)/pi;
 
  ay=1/20;
 xamax=30;
 qqy=(xamax)^2*(4-pi)/pi;
 
xe=[0 0 0 0 0 0]';
 p=100000*eye(6);
 
M=10;
 
ap=4;v=4;
R=(0.2303*ap/v)^2;
 
%%%估计横轴、纵轴使用最小二乘估计EKF方法， 
%[xx1,xxe1,P33,NN]=funDataDrivenModelYWwithEKF(TT,R,ax,qqx,ay,qqy,xe,p,dm,N,readerxy);
 
%%%%使用UKF计算的方法，可以节约计算时间，估计横轴、纵轴使用最小二乘估计Yule-Walker方法
[xx1,xxe1,P33,NN,qqxx,RR0x,RR0y]=funDataDrivenModelYWnonparell(TT,R,ax,qqx,ay,qqy,xe,p,dm,N,readerxy);

cov1=xys-[xx1(1,:);xx1(4,:)];%计算方差
mm=mean(cov1,2);
cov1=[cov1(1,:)-mm(1);cov1(2,:)-mm(2)];
cov1=diag(cov1*cov1');
cov1=cov1/length(xys(1,:));
covvv=[covvv cov1];
 
end
 
covvv;
XY=sum(covvv,2)/N   
 %%%%%画出结果
XY2=sqrt(XY(1)*XY(1)+XY(2)*XY(2))
plot(xys(1,:),xys(2,:),'k-.','LineWidth',4);hold on
plot(xx1(1,:),xx1(4,:),'r*');hold off
legend('the reference trajectory','the estimation trajectory')
  figure 
 
 subplot(2,1,1),plot(ts,xys(1,:),'k-.','LineWidth',4);hold on,plot(ts,xx1(1,:),'r*')
 legend('the reference trajectory','the estimation trajectory')
 xlabel('time'),ylabel('Horizontal axis tracking')
 
subplot(2,1,2),plot(ts,xys(2,:),'k-.','LineWidth',4);hold on,plot(ts,xx1(4,:),'r*')
legend('the reference trajectory','the estimation trajectory')
xlabel('time'),ylabel('Longitudinal axis tracking')   
 
 figure
 
subplot(2,1,1),plot(ts,xys(1,:)-xx1(1,:),'k','LineWidth',3)
xlabel('time'),ylabel('Error of horizontal axis tracking')
subplot(2,1,2),plot(ts,xys(2,:)-xx1(4,:),'k','LineWidth',3)
xlabel('time'),ylabel('Error of longitudinal axis tracking')  
