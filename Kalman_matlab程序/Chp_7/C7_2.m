clc
clear
%%%�������������ںͲ������������
%ʹ��kalman�˲��������λ������й���
%%%�������λ����Ĳ���ֵ
T=0.3;
t=1:T:30;a=10;omig=pi/3;R=2;
[yreal,ym]=funtrackingsnake(a,omig,t,R);
%%%����
qq=10;%%%��ֵ��ѡ�������Ӱ��ܴ󣬿���ȥ��ͬ��ֵ����ȡ0.1��1��10�ȣ�ֱ���޸ĳ����е�ֵ���ɡ�

[A,Q]=CAmodel(T,qq);C=[1 0 0];

xe=zeros(length(Q),1);p=1000*eye(size(A));xx1=[];
for i=1:length(t)
[xe,p]=kalmanfun(A,C,Q,R,xe,ym(i),p)
xx1=[xx1 xe];
end
plot(t,ym,t,C*xx1)

%%%%%%%%%%ѡ��ģ��
