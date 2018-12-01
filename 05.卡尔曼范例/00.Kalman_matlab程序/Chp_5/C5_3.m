clc
clear
b=0.01;d=0.35;
x(:,1)=[1+sqrt(2)*b;pi/2];
x(:,2)=[1;pi/2+sqrt(2)*d];
x(:,3)=[1-sqrt(2)*b;pi/2];
x(:,4)=[1;pi/2-sqrt(2)*d];
 
 %下面利用循环完成（5.42b），求非线性.
for i=1:4
    z(:,i)=[x(1,i)*cos(x(2,i));x(1,i)*sin(x(2,i))];
end
%我们先把所求的非线性变量画出来.
hold on
plot(z(1,1),z(2,1),'ko','MarkerSize',8,'MarkerFaceColor',[.5 1 .3]);
plot(z(1,2),z(2,2),'ko','MarkerSize',8,'MarkerFaceColor',[.5 1 .3]);
plot(z(1,3),z(2,3),'ko','MarkerSize',8,'MarkerFaceColor',[.5 1 .3]);
plot(z(1,4),z(2,4),'ko','MarkerSize',8,'MarkerFaceColor',[.5 1 .3]);

%根据（5.43）式求平均值，再把均值画出来.
sigmamean=sum(z,2)/4;
plot(sigmamean(1),sigmamean(2),'^','MarkerSize',8,'MarkerFaceColor',[.5 1 .3])
 %根据（5.44）求方差.
Pz=zeros(2);
 for i=1:4
Pz=Pz+(z(:,i)-sigmamean)*(z(:,i)-sigmamean)';
 end
%将求得的方差值求平方后，在画出一个范围，是一个椭圆.
Pz=sqrt(Pz/4);
rectangle('Position',[sigmamean(1)-Pz(1,1),sigmamean(2)-Pz(2,2),Pz(1,1)*2,Pz(2,2)*2],'Curvature',[1,1],'LineWidth',3)