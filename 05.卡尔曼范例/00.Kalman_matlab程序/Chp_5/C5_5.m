clc
clear
n=2;nn=4;
W0=0;
W(1)=(1-W0)/(n+1);
W(2)=W(1);
W(3)=W(1);
 
gama(:,1)=[0;0];
gama(:,2)=[-1/sqrt(2*W(1));-1/sqrt(6*W(1))];
gama(:,3)=[1/sqrt(2*W(1));-1/sqrt(6*W(1))];
gama(:,4)=[0;2/sqrt(6*W(1))];
    
meanx=[1;pi/2];
Px=diag([0.01^2,0.35^2]);
for i=1:4
x(:,i)=meanx+sqrt(Px)*gama(:,i);
end
 
for i=1:nn
    y(:,i)=[x(1,i)*cos(x(2,i));x(1,i)*sin(x(2,i))];
end
 
sigmamean=W0*y(:,1);
for i=1:nn-1
sigmamean=sigmamean+W(i)*y(:,i+1);
end
 
hold on
plot(y(1,1),y(2,1),'ko','MarkerSize',8,'MarkerFaceColor',[.5 0.90 .9]);
plot(y(1,2),y(2,2),'ko','MarkerSize',8,'MarkerFaceColor',[.5 0.90 .9]);
plot(y(1,3),y(2,3),'ko','MarkerSize',8,'MarkerFaceColor',[.5 0.90 .9]);
plot(y(1,4),y(2,4),'ko','MarkerSize',8,'MarkerFaceColor',[.5 0.90 .9]);
 
plot(sigmamean(1),sigmamean(2),'^','MarkerSize',8,'MarkerFaceColor',[.5 1 .3])
 
Py=W0*(y(:,1)-sigmamean)*(y(:,1)-sigmamean)';
for i=1:3
Py=Py+W(i)*(y(:,i+1)-sigmamean)*(y(:,i+1)-sigmamean)';
end
Py=sqrt(Py);
rectangle('Position',[sigmamean(1)-Py(1,1),sigmamean(2)-Py(2,2),Py(1,1)*2,Py(2,2)*2],'Curvature',[1,1],'LineWidth',1)