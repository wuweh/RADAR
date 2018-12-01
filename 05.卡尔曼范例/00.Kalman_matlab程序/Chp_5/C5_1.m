clc
clear
a=-0.01;b=0.01;
rnoise= a + (b-a).*rand(300,1);
rmean=1;
r=rnoise+rmean; 
c=-0.35;d=0.35;
sitanoise= c + (d-c).*rand(300,1);
sitamean=pi/2;
sita=sitanoise+sitamean;
y1=r.*cos(sita);
y2=r.*sin(sita);
meany1=mean(y1);
meany2=mean(y2);
%»­³ö½á¹ûÍ¼.
plot(y1,y2,'.')
xlabel('y1');ylabel('y2')
hold on
plot(0,1,'ko')
plot(meany1,meany2,'k*')

