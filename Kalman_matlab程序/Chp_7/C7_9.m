clc
clear
summ=0;N=10;covvv=[];
for n=1:N
 
%%?��ģ��켣target���и���
 
 load mytarget 
%%%%���ò������ڡ��������̵Ĳ���
T=0.1;
R=25; 
C=[1 0 0];

y=xys+sqrt(R)*randn(size(xys));
%%%%%%%%%%%���ƺ��ᡢ����ģ�Ͳ�����ֵ
 a=1/20;
 xamax=5;
 qq=(xamax)^2*(4-pi)/pi;
 [RR0,RR1,xx1,xxe1,aa1,qq1,Ea1]=StartrackingModel(y(1,:),a,qq,T,R,C);
%%%%�ڹ��ƺ����ʱ��ģ�Ͳ����ĳ�ֵ�����˱仯����ˣ��ڹ¼������ʱ���ٴ��������ֵ
a=1/20;
 xamax=5;
 qq=(xamax)^2*(4-pi)/pi;
[RR02,RR12,xx2,xxe2,aa2,qq2,Ea2]=StartrackingModel(y(2,:),a,qq,T,R,C);
 
%%%%���ƽ����󣬼��㷽��ͻ����ͼ��
covv=diag(cov(xys'-[C*xx1;C*xx2]'))
covvv=[covvv covv];
summ=summ+covv;
end
summ/N
 
plot(xys(1,:),xys(2,:),'b-.');hold on
plot(C*xx1,C*xx2,'--');hold off
legend('the real trajectory','the estimation trajectory')
 figure 
 subplot(2,1,1),plot(ts,xys(1,:),'-',ts,C*xx1,'--')
 legend('the real trajectory','the estimation trajectory')
 xlabel('time'),ylabel('Horizontal axis tracking')
 subplot(2,1,2),plot(ts,xys(2,:),'-',ts,C*xx2,'--')
legend('the real trajectory','the estimation trajectory')
xlabel('time'),ylabel('Longitudinal axis tracking')   
figure
subplot(2,1,1),plot(ts,xys(1,:)-C*xx1)
xlabel('time'),ylabel('Horizontal axis tracking')
subplot(2,1,2),plot(ts,xys(2,:)-C*xx2)
xlabel('time'),ylabel('Longitudinal axis tracking')   
figure
subplot(3,1,1),plot(ts(2:length(ts)),aa1)
xlabel('(a)')
subplot(3,1,2),plot(ts(2:length(ts)),qq1)
xlabel('(b)')
subplot(3,1,3),plot(ts(2:length(ts)),Ea1)
xlabel('(c)')
figure
subplot(3,1,1),plot(ts(2:length(ts)),aa2)
xlabel('(a)')
subplot(3,1,2),plot(ts(2:length(ts)),qq2)
xlabel('(b)')
subplot(3,1,3),plot(ts(2:length(ts)),Ea2)
xlabel('(c)')