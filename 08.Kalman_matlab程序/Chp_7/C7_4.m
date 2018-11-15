clc
clear
number=100;
x=zeros(1,100);y=zeros(1,number);
ax1=0.075;ay1=0.075;
ax2=0.3;ay2=-0.3;
T=10;
t=0:T:T*(number-1);
for i=1:number
    if i==1
        x(i)=2000;
           y(i)=10000;
    elseif i<=39
        x(i)=x(i-1);
         y(i)=y(i-1)-5;
    elseif i<=60;
        x(i)=2*x(i-1)-x(i-2)+ax1;
         y(i)=2*y(i-1)-y(i-2)+ay1;
    elseif i<=66;
         x(i)=2*x(i-1)-x(i-2)+ax2;
          y(i)=2*y(i-1)-y(i-2)+ay2;
    else
         x(i)=x(i-1);
          y(i)=y(i-1)-5;
    end
    
end
ts=t;
xys=[x;y];
subplot(2,1,1);plot(t,x)
subplot(2,1,2);plot(t,y)
figure
plot(x,y)