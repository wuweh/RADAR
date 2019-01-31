clc
clear all
close all

tic
T = 1;
n=50;
global Pd Pg C lambda A Q R ;
Pd=1;       
Pg=0.9997;   

ganmma = 16;
g_sigma = 9.21;     
lambda =0.0004;          

target_position=zeros(4,n);     

x_filter=zeros(4,n); 
x_filter1=zeros(4,n);

A = [1 T 0 0;
     0 1 0 0;
     0 0 1 T;
     0 0 0 1];  
 
Q = [ 0.01 0 0 0; 
      0 0.01 0 0; 
      0 0 0.01 0; 
      0 0 0 0.01];     
 
C = [1 0 0 0;
     0 0 1 0];     

target_delta=1;           
R=[ 1^2 0;
    0 1^2]; 

R11=target_delta; 
R22=target_delta; 
R12=0; 
R21=0;

P=[R11 R11/T R12 R12/T;     R11/T 2*R11/T^2 R12/T 2*R12/T^2; 
   R21 R21/T R22 R22/T;     R21/T 2*R21/T^2 R22/T 2*R22/T^2];    

X0=[200;0;500;15];            
target_position(:,1)=X0; 
Vk=[target_delta*randn;
    target_delta*randn]; 
Zk(:,1)=C*target_position(:,1)+Vk;   


for i=2:1:n 
    target_position(:,i)=A*target_position(:,i-1); 
    Vk=[target_delta*randn;
        target_delta*randn]; 
    Zk(:,i)=C*target_position(:,i)+Vk;     
end

Noise=[];
NOISE=[]; 

x_filter(:,1)  = X0;
sumx = 0;
figure;
for t=2:n
    x_predic = A*x_filter(:,t-1);  

    [P_predic, Z_predic, S, K] = kalman_part_func(x_predic, P);
    
    Av=pi*ganmma*sqrt(det(S))   ;                         
    number_returns=floor(10*Av*lambda+1);     %计算杂波个数          
    side=sqrt(10*Av)/2;                                    
    Noise_x=target_position(1,t)-side+2*rand(1,number_returns)*side;            
    Noise_y=target_position(3,t)-side+2*rand(1,number_returns)*side; 
    Noise=[Noise_x ;Noise_y]; 
     
    b=zeros(1,2); 
    b(1)=Zk(1,t); 
    b(2)=Zk(2,t); 

    y1=[Noise b'];   
    y=[]; 
    d=[]; 
    m=0; 
    
    %判断量测点是否与当前航迹相关（卡方值）
    for j=1:(number_returns+1)     
        d=y1(:,j)-Z_predic; 
        D=d'*inv(S)*d;  
        if D<=16 
            plot(y1(1,j),y1(2,j),'.');
            hold on;
            y=[y y1(:,j)];   
            m=m+1;          
        end 
    end 
    m
    %% PDA_Function
    [x_putput, P] = PDA_Function(x_predic, P_predic, S, Z_predic, K,y, m);
    x_filter(:,t) = x_putput;
    
end

plot(target_position(1,:),target_position(3,:),'-bo');hold on; 
plot(Zk(1,:),Zk(2,:),'g*') ;hold on; 
plot(x_filter(1,:),x_filter(3,:),'-r*');hold on; 
axis([190 210 400 1400]);
toc

