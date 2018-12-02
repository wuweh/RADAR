clc
clear all
close all
summ=0;
N=1;
for n=1:N
    load mytarget 
    T=0.1;
    R=1;
    
    y=xys+sqrt(R)*randn(size(xys));
    pa=[];

    %%%IMM模型中使用的3个基本模型  
    %匀加速模型（CA_Model）
    [A1,Q1]=CAmodel(T,0.1);
    [A2,Q2]=CAmodel(T,1);
    [A3,Q3]=CAmodel(T,10);
    C=[1 0 0];

    %%%模型转移参数
    Pij=[0.95 0.05 0;0.33 0.34 0.33;0 0.05 0.95];
    model0=[1/3,1/3,1/3];

    %%%%估计横轴
    %%%%%%%%%设置状态初值
    x0=[0 0 0]';
    [xx1]=IMM(A1,Q1,A2,Q2,A3,Q3,C,R,Pij,model0,y(1,:),x0);

    % %%%%估计纵轴
    [xx2]=IMM(A1,Q1,A2,Q2,A3,Q3,C,R,Pij,model0,y(2,:),x0);
end
 
plot(xys(1,:),xys(2,:),'b-');hold on
plot(y(1,:),y(2,:),'.');hold on
plot(C*xx1,C*xx2,'r-.');hold on

figure 
subplot(2,1,1),plot(ts,xys(1,:),ts,C*xx1)
subplot(2,1,2),plot(ts,xys(2,:),ts,C*xx2)
figure 
subplot(2,1,1),plot(ts,xys(1,:)-C*xx1)
subplot(2,1,2),plot(ts,xys(2,:)-C*xx2)
covv=diag(cov(xys'-[C*xx1;C*xx2]'));
% 
summ=summ+covv;
summ/N;

