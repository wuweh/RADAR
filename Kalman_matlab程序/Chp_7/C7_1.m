clc
clear
%%%�������������ںͲ������������
T=1;
R=40^2;;
%%%Բ�λ����Ĳ���ֵ
Tt=31;Tz=57;[t,y]=funtrackinglinecircle(T,Tt,Tz,R);
%%%%%%%%%%ѡ��ģ��
qq=1;%%%��ֵ��ѡ�������Ӱ��ܴ󣬿���ȥ��ͬ��ֵ����ȡ0.1��1��10�ȣ�ֱ���޸ĳ����е�ֵ���ɡ�

 [A,Q]=CVmodel(T,qq);C=[1 0]; 
% [A,Q]=CAmodel(T,qq);C=[1 0 0]; %%������ģ����ִ�е�ʱ��ѡ���е�һ����
%ʹ��Kalman�����˲�,����Բ�λ�������Ҫ������ֱ����
%%%���ƺ���
xe=zeros(length(Q),1);p=1000*eye(size(A));xx1=[];
for i=1:length(t)
[xe,p]=kalmanfun(A,C,Q,R,xe,y(1,i),p);
xx1=[xx1 xe];
end
%%%%��������
xe=zeros(length(Q),1);p=1000*eye(size(A));xx2=[];
for i=1:length(t)
[xe,p]=kalmanfun(A,C,Q,R,xe,y(2,i),p)
xx2=[xx2 xe];
end
plot(y(1,:),y(2,:),'*');hold on,plot(C*xx1,C*xx2);hold off
figure