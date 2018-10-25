clc
clear
t=0.01:0.01:10;
za=[20+sqrt(2)*randn(500,1);30+sqrt(2)*randn(500,1)];
h=1;
w=0.5;
M=1000;                     
ea0=0;
for i=1:1000
M=inv(inv(M)+h'*w*h);
ea(i)=ea0+M*h'*w*(za(i)-h*ea0);
ea0=ea(i);
end
plot(ea)