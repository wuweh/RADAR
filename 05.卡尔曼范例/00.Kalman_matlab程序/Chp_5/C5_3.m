clc
clear
b=0.01;d=0.35;
x(:,1)=[1+sqrt(2)*b;pi/2];
x(:,2)=[1;pi/2+sqrt(2)*d];
x(:,3)=[1-sqrt(2)*b;pi/2];
x(:,4)=[1;pi/2-sqrt(2)*d];
 
 %��������ѭ����ɣ�5.42b�����������.
for i=1:4
    z(:,i)=[x(1,i)*cos(x(2,i));x(1,i)*sin(x(2,i))];
end
%�����Ȱ�����ķ����Ա���������.
hold on
plot(z(1,1),z(2,1),'ko','MarkerSize',8,'MarkerFaceColor',[.5 1 .3]);
plot(z(1,2),z(2,2),'ko','MarkerSize',8,'MarkerFaceColor',[.5 1 .3]);
plot(z(1,3),z(2,3),'ko','MarkerSize',8,'MarkerFaceColor',[.5 1 .3]);
plot(z(1,4),z(2,4),'ko','MarkerSize',8,'MarkerFaceColor',[.5 1 .3]);

%���ݣ�5.43��ʽ��ƽ��ֵ���ٰѾ�ֵ������.
sigmamean=sum(z,2)/4;
plot(sigmamean(1),sigmamean(2),'^','MarkerSize',8,'MarkerFaceColor',[.5 1 .3])
 %���ݣ�5.44���󷽲�.
Pz=zeros(2);
 for i=1:4
Pz=Pz+(z(:,i)-sigmamean)*(z(:,i)-sigmamean)';
 end
%����õķ���ֵ��ƽ�����ڻ���һ����Χ����һ����Բ.
Pz=sqrt(Pz/4);
rectangle('Position',[sigmamean(1)-Pz(1,1),sigmamean(2)-Pz(2,2),Pz(1,1)*2,Pz(2,2)*2],'Curvature',[1,1],'LineWidth',3)