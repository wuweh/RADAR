clc
clear
summ=0;N=10;
for n=1:N

%%读入mytarget 数据文件
load mytarget 
T=0.1;
R=25;
 
y=xys+sqrt(R)*randn(size(xys));
qqa=[];
%%%%%%%%%%%选择模型参数
 a=1/20;
 xamax=30;
 C=[1 0 0];
%%%估计横轴
xe=zeros(3,1);p=10*eye(3);xx1=[];
for i=1:length(y(1,:))
xa=xe(3);
[A1,A,Q,U,qa]=Starmodel(T,xa,a,xamax);
[xe,p]=kalmanadfun(A1,A,U,C,Q,R,xe,y(1,i),p);
xa=xe(3);
xx1=[xx1 xe];
qqa=[qqa qa];
end
%%%%估计纵轴
xe=zeros(3,1);p=10*eye(3);xx2=[];
for i=1:length(y(2,:))
xa=xe(3);
[A1,A,Q,U,qa]=Starmodel(T,xa,a,xamax);
[xe,p]=kalmanadfun(A1,A,U,C,Q,R,xe,y(2,i),p);
xx2=[xx2 xe];
end
 
covv=diag(cov(xys'-[C*xx1;C*xx2]'))

summ=summ+covv;
end
summ/N
 plot(xys(1,:),xys(2,:),'--');hold on
plot(C*xx1,C*xx2,'r*');hold off
figure 
subplot(2,1,1),plot(ts,xys(1,:),ts,C*xx1,'-.')
subplot(2,1,2),plot(ts,xys(2,:),ts,C*xx2,'-.')
figure 
subplot(2,1,1),plot(ts,xys(1,:)-C*xx1)
subplot(2,1,2),plot(ts,xys(2,:)-C*xx2)