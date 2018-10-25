clc
clear
%%%%%%%%%%%%%%%%%%�������������� p(k+1)=a*p(k)+w(t)
a=0.8;   
q=2;          %%w(t)�ķ���
p0=3;p=[];
for i=1:4000
    p=[p p0];
    p0=a*p0+sqrt(q)*randn(1);
end
 
sum_ar=0;
sum_qr=0;
sum_aY=0;
sum_qY=0;

J=500;
for j=1:J  %%%%%%%%%%%%����500�ι���������ƽ�����
%%%%%%%%%%%%%%%%%%%%%%%%������С���˵��Ʒ����õ��Ĺ��Ʋ���
[ar,qr]=LSRecursive(p);
arlast=ar(max(length(ar)));
qrlast=qr(max(length(qr)));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����Yule-Walker�����õ��Ĺ��Ʋ���
[aY,qY]=YuleWalker(p);
aYlast=aY(max(length(aY)));
qYlast=qY(max(length(qY)));
 
sum_ar=sum_ar+arlast;
sum_qr=sum_qr+qrlast;
sum_aY=sum_aY+aYlast;
sum_qY=sum_qY+qYlast;
end
 
arl=sum_ar/J%%%%��ƽ��ֵ
qrl=sum_qr/J
aYl=sum_aY/J
qYl=sum_qY/J
 
subplot(2,1,1),plot(ar,'k','LineWidth',2)
xlabel('k'),ylabel('����ϵ��a�Ĺ���')
subplot(2,1,2),plot(qr,'k','LineWidth',2)
xlabel('k'),ylabel('�������')
figure
subplot(2,1,1),plot(aY,'k','LineWidth',2)
xlabel('k'),ylabel('����ϵ��a�Ĺ���')
subplot(2,1,2),plot(qY,'k','LineWidth',2)
xlabel('k'),ylabel('�������')