clc
clear
summ=0;N=10;
for n=1:N

%%����Hotarget�����ļ�
load Hotarget 
T=10;
R=10000 %ts xys;
 
y=xys+sqrt(R)*randn(size(xys));
qqa=[];
%%%%%%%%%%%ѡ��ģ�Ͳ���
 a=1/20;
 xamax=30;
 C=[1 0 0];
%%%���ƺ���
xe=zeros(3,1);p=10*eye(3);xx1=[];
for i=1:length(y(1,:))
xa=xe(3);
[A1,A,Q,U,qa]=Starmodel(T,xa,a,xamax);
[xe,p]=kalmanadfun(A1,A,U,C,Q,R,xe,y(1,i),p);
xa=xe(3);
xx1=[xx1 xe];
qqa=[qqa qa];
end
%%%%��������
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