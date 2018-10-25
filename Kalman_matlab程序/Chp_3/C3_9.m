clc;
clear;
t=0:0.01:9.99;
m=2*t+20;
v=2;
s3=m+sqrt(v)*randn(1,1000);
plot(t,s3)
xlabel('t')
ylabel('measurement data')
ea=[0;0];           
M=diag(1)*1000;                   
estiamte=[];
W=0.5;
 
for i=1:1000
h=[t(i) 1];
M=inv(inv(M)+h'*W*h);
ea=ea+M*h'*W*(s3(i)-h*ea);
 
estiamte=[estiamte,ea];
end
 
plot(estiamte(1,:));figure
plot(estiamte(2,:));