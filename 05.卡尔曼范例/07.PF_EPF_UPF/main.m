%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 主函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 功能说明：ekf,ukf,pf,epf,upf算法的综合比较程序
function main
clc;close all;clear all;
rand('seed',3);
randn('seed',6);
%error('下面的参数T请参考书中的值设置，然后删除本行代码') 
T = 50;
R =  1e-5;              
                       
g1 = 3;                
g2 = 2;               
 
X = zeros(1,T);
Z = zeros(1,T);
processNoise = zeros(T,1);
measureNoise = zeros(T,1);
X(1) = 1;                         
P0 = 3/4;
 
Qekf=10*3/4;                
Rekf=1e-1;                    
Xekf=zeros(1,T);               
Pekf = P0*ones(1,T);         
Tekf=zeros(1,T);               

 
Qukf=2*3/4;                  
Rukf=1e-1;                    
Xukf=zeros(1,T);               
Pukf = P0*ones(1,T);          
Tukf=zeros(1,T);                

 
N=200;                      
Xpf=zeros(1,T);             
Xpfset=ones(T,N);         
Tpf=zeros(1,T);            

 
Xepf=zeros(1,T);           
Xepfset=ones(T,N);         
Pepf = P0*ones(T,N);          
Tepf=zeros(1,T);               

 
Xupf=zeros(1,T);            
Xupfset=ones(T,N);              
Pupf = P0*ones(T,N);       
Tupf=zeros(1,T);               
 
for t=2:T
    processNoise(t) =  gengamma(g1,g2);  
    measureNoise(t) =  sqrt(R)*randn;    
 
    X(t) = feval('ffun',X(t-1),t) +processNoise(t);
    Z(t) = feval('hfun',X(t),t) + measureNoise(t);
    
    tic
    [Xekf(t),Pekf(t)]=ekf(Xekf(t-1),Z(t),Pekf(t-1),t,Qekf,Rekf);
    Tekf(t)=toc;
    
 
    tic
    [Xukf(t),Pukf(t)]=ukf(Xukf(t-1),Z(t),Pukf(t-1),Qukf,Rukf,t);
    Tukf(t)=toc;
    
 
    tic
    [Xpf(t),Xpfset(t,:)]=pf(Xpfset(t-1,:),Z(t),N,t,R,g1,g2);
    Tpf(t)=toc;
    
 
    tic
    [Xepf(t),Xepfset(t,:),Pepf(t,:)]=epf(Xepfset(t-1,:),Z(t),t,Pepf(t-1,:),N,R,Qekf,Rekf,g1,g2);
    Tepf(t)=toc;
    
 
    tic
    [Xupf(t),Xupfset(t,:),Pupf(t,:)]=upf(Xupfset(t-1,:),Z(t),t,Pupf(t-1,:),N,R,Qukf,Rukf,g1,g2);
    Tupf(t)=toc;
end
 
ErrorEkf=abs(Xekf-X); 
ErrorUkf=abs(Xukf-X);  
ErrorPf=abs(Xpf-X);     
ErrorEpf=abs(Xepf-X);   
ErrorUpf=abs(Xupf-X);   
 
figure
hold on;box on;
p1=plot(1:T,X,'-k.','lineWidth',2);
p2=plot(1:T,Xekf,'m:','lineWidth',2);
p3=plot(1:T,Xukf,'--','lineWidth',2);
p4=plot(1:T,Xpf,'-ro','lineWidth',2);
p5=plot(1:T,Xepf,'-g*','lineWidth',2);
p6=plot(1:T,Xupf,'-b^','lineWidth',2);
legend([p1,p2,p3,p4,p5,p6],'真实状态','EKF估计','UKF估计','PF估计','EPF估计','UPF估计')
xlabel('Time','fontsize',10)
title('Filter estimates (posterior means) vs. True state','fontsize',10)

 
figure
hold on;box on;
p1=plot(1:T,ErrorEkf,'-k.','lineWidth',2);
p2=plot(1:T,ErrorUkf,'-m^','lineWidth',2);
p3=plot(1:T,ErrorPf,'-ro','lineWidth',2);
p4=plot(1:T,ErrorEpf,'-g*','lineWidth',2);
p5=plot(1:T,ErrorUpf,'-bd','lineWidth',2);
legend([p1,p2,p3,p4,p5],'EKF偏差','UKF偏差','PF偏差','EPF偏差','UPF偏差')

 
figure
hold on;box on;
p1=plot(1:T,Tekf,'-k.','lineWidth',2);
p2=plot(1:T,Tukf,'-m^','lineWidth',2);
p3=plot(1:T,Tpf,'-ro','lineWidth',2);
p4=plot(1:T,Tepf,'-g*','lineWidth',2);
p5=plot(1:T,Tupf,'-bd','lineWidth',2);
legend([p1,p2,p3,p4,p5],'EKF时间','UKF时间','PF时间','EPF时间','UPF时间')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%