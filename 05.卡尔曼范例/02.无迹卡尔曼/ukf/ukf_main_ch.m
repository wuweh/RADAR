%n=3; %number of state
clc;
clear all;
close all;
n=6;
t=0.2;
q=0.3; %std of process
r=0.5; %std of measurement
Q=q^2*eye(n); % covariance of process 
R=r^1; % covariance of measurement 
%f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))]; % nonlinear state equations
%实际方程
f=@(x)[x(1)+t*x(3);x(2)+t*x(4);x(3)+t*x(5);x(4)+t*x(6);x(5);x(6)]; % nonlinear state equations
h=@(x)[sqrt(x(1)+1);0.8*x(2)+0.3*x(1);x(3);x(4);x(5);x(6)]; %
% measurement equation
%s=[0;0;1]; % initial state
s=[0.3;0.2;1;2;2;-1];
x=s+q*randn(n,1); %initial state with noise 
P = eye(n); % initial state covraiance
N=100; % total dynamic steps
xV = zeros(n,N); %estmate % allocate memory
sV = zeros(n,N); %actual
zV = zeros(n,N);
for k=1:N
    z = h(s) + r*randn; % measurments  
    sV(:,k)= s; % save actual state  
    zV(:,k) = z; % save measurment 
    [x, P] = ukf(f,x,P,h,z,Q,R); % ukf
    xV(:,k) = x; % save estimate %
    s = f(s) + q*randn(n,1); % update process 
end
for k=1:4 % plot results
    subplot(4,1,k)
    plot(1:N, sV(k,:), '-r', 1:N, xV(k,:), '--b',1:N,zV(k,:),'*g')
end