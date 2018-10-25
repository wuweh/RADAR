clc
clear
b=0.01;d=0.35;

k=1;n=2;
x(:,1)=[1+sqrt(k+n)*b;pi/2];
x(:,2)=[1;pi/2+sqrt(k+n)*d];
x(:,3)=[1-sqrt(k+n)*b;pi/2];
x(:,4)=[1;pi/2-sqrt(k+n)*d];
x0=[1,pi/2];
 
nn=4;
w(1:nn)=1/(2*(n+k));
w0=k/(k+n);
 
y0=[x0(1)*cos(x0(2));x0(1)*sin(x0(2))];
for i=1:nn
    y(:,i)=[x(1,i)*cos(x(2,i));x(1,i)*sin(x(2,i))];
end
 
sigmamean=w0*y0;
for i=1:nn
sigmamean=sigmamean+w(i)*y(:,i);
end
 
hold on

plot(y0(1),y0(2),'ko','MarkerSize',8,'MarkerFaceColor',[.5 0.10 .9]);
plot(y(1,1),y(2,1),'ko','MarkerSize',8,'MarkerFaceColor',[.5 0.10 .9]);
plot(y(1,2),y(2,2),'ko','MarkerSize',8,'MarkerFaceColor',[.5 0.10 .9]);
plot(y(1,3),y(2,3),'ko','MarkerSize',8,'MarkerFaceColor',[.5 0.10 .9]);
plot(y(1,4),y(2,4),'ko','MarkerSize',8,'MarkerFaceColor',[.5 0.10 .9]);
 
plot(sigmamean(1),sigmamean(2),'^','MarkerSize',8,'MarkerFaceColor',[.5 1 .3])

Py=w0*(y0-sigmamean)*(y0-sigmamean)';
for i=1:4
Py=Py+w(i)*(y(:,i)-sigmamean)*(y(:,i)-sigmamean)';
end
Py=sqrt(Py);
rectangle('Position',[sigmamean(1)-Py(1,1),sigmamean(2)-Py(2,2),Py(1,1)*2,Py(2,2)*2],'Curvature',[1,1],'LineWidth',3)