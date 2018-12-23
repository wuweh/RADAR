clear all;
clc;
L = 10;
T = 1;
NT = 100;               %Á£×ÓÊýÄ¿
v = 0.10;
r1 = 16;
r2 = 0.002^2;
R = diag([r1,r2]);
Q1 =diag([1,0.01,1,0.01]); 
A1 = [1,T,0,0;
      0,1,0,0;
      0,0,1,T;
      0,0,0,1];
  
Q2 = diag([2,0.01,1,0.01]); 
A2 = [1,sin(v*T)/v,0,-(1-cos(v*T))/v;
     0,cos(v*T),0 -sin(v*T);
     0,(1-cos(v*T))/v,1,sin(v*T)/v;
     0,sin(v*T),0,cos(v*T)];
 
x(:,1)  = [200,15,150,6]';
xe(:,1) = [150,10,100,5]';

for t = 2:34
     x(:,t) = A1*x(:,t-1)+sqrtm(Q1)*randn(4,1);
end
for t = 35:67
     x(:,t) = A2*x(:,t-1)+sqrtm(Q2)*randn(4,1);
end
for t = 68:100
     x(:,t) = A1*x(:,t-1)+sqrtm(Q1)*randn(4,1);
end
for t = 1:NT
     r(t) = sqrt(x(1,t)^2+x(3,t)^2)+sqrt(r1)*randn;
     sita(t) = atan(x(3,t)./x(1,t))+sqrt(r2)*randn; 
     y(:,t) = [r(t);sita(t)];
              
end



E1 = zeros(1,NT);
E1x = zeros(1,NT);
E1y = zeros(1,NT);

E2 = zeros(1,NT);
E2x = zeros(1,NT);
E2y = zeros(1,NT);

E3 = zeros(1,NT);
E3x = zeros(1,NT);
E3y = zeros(1,NT);

for n = 1:L
    [X1] = immpf1(x,y,R,xe(:,1));
    [X2] = immpf2(x,y,R,xe(:,1));
    [X3] = immpf3(x,y,R,xe(:,1));
    
    for k = 1:NT
         E1(k)  = E1(k)+((X1(1,k)-x(1,k))^2+(X1(3,k)-x(3,k))^2);
         E1x(k) = E1x(k)+(X1(1,k)-x(1,k))^2;
         E1y(k) = E1y(k)+(X1(3,k)-x(3,k))^2;
         
         E2(k)  = E2(k)+((X2(1,k)-x(1,k))^2+(X2(3,k)-x(3,k))^2);
         E2x(k) = E2x(k)+(X2(1,k)-x(1,k))^2;
         E2y(k) = E2y(k)+(X2(3,k)-x(3,k))^2;
         
         E3(k)  = E3(k)+((X3(1,k)-x(1,k))^2+(X3(3,k)-x(3,k))^2);
         E3x(k) = E3x(k)+(X3(1,k)-x(1,k))^2;
         E3y(k) = E3y(k)+(X3(3,k)-x(3,k))^2;
    end
    n
end

figure(1)
plot(x(1,:),x(3,:),'k');hold on;grid on;
plot(X1(1,:),X1(3,:),'b');
plot(X2(1,:),X2(3,:),'g');
plot(X3(1,:),X3(3,:),'r');


E1 = sqrt(E1/L);
E2 = sqrt(E2/L);
E3 = sqrt(E3/L);

E1x = sqrt(E1x/L);
E2x = sqrt(E2x/L);
E3x = sqrt(E3x/L);

E1y = sqrt(E1y/L);
E2y = sqrt(E2y/L);
E3y = sqrt(E3y/L);

t = 1:NT;
figure(2)
plot(t,E1,'k:',t,E2,'k--',t,E3,'k');
ylabel('Distance error(m)');
xlabel('Time(s)');
legend('IMMPF1', 'IMMPF2','IMMPF3'); 

figure(3)
plot(t,E1x,'r',t,E2x,'b',t,E3x,'k');
ylabel('X error(m)');
xlabel('Time(s)');
legend('IMMPF1', 'IMMPF2','IMMPF3'); 

figure(4)
plot(t,E1y,'r',t,E2y,'b',t,E3y,'k');
ylabel('Y error(m)');
xlabel('Time(s)');
legend('IMMPF1', 'IMMPF2','IMMPF3'); 



RMS1 = sum(E1)/NT;
RMS2 = sum(E2)/NT;
RMS3 = sum(E3)/NT;