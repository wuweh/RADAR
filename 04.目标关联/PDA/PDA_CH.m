clc
clear all
close all

tic
T = 1;
n=50;
global Pd Pg gamma C;
Pd=1;       
Pg=0.99;   

g_sigma = 9.21;     
gamma = 1;          

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
     0 1 0 0
     0 0 1 0;
     0 0 0 1];     

target_delta=1;           
R=[ 1^2 0 0 0;
    0 0.1^2 0 0;
    0 0 1^2 0 ;
    0 0 0 0.1^2]; 

% G=[T^2/2 0;
%    0 T^2/2;
%    T 0;
%    0 T]; 

R11=target_delta; 
R22=target_delta; 
R12=0; 
R21=0;

P=[R11 R11/T R12 R12/T;     R11/T 2*R11/T^2 R12/T 2*R12/T^2; 
   R21 R21/T R22 R22/T;     R21/T 2*R21/T^2 R22/T 2*R22/T^2];    

X0=[200;0;500;15];            
target_position(:,1)=X0; 
Vk=[target_delta*randn;
    0.1^2;
    target_delta*randn;
    0.1^2]; 
Zk(:,1)=C*target_position(:,1)+Vk;   


for i=2:1:n 
    target_position(:,i)=A*target_position(:,i-1); 
    
%     Vk=[target_delta*randn;target_delta*randn];  
    Vk=[target_delta*randn;
        0.1^2;
        target_delta*randn;
        0.1^2]; 
    Zk(:,i)=C*target_position(:,i)+Vk;     
end

Noise=[];
NOISE=[]; 

sumx = 0;
figure;
for t=1:n
    if t~=1 
        x_predic = A*x_filter(:,t-1);  
    else 
        x_predic = target_position(:,1);   
    end 
    
    plot(Zk(1,t),Zk(3,t),'g*') ;
    hold on; 

    P_predic = A*P*A'+Q; 
    Z_predic = C*x_predic;
    
    S = C*P_predic*C'+ R; 
    K = P_predic*C'*inv(S);
    
    Av=pi*g_sigma*sqrt(det(S))   ;                         
%     number_returns=floor(10*Av*gamma+1) ;                
    number_returns = 20;
    side=sqrt(10*Av)/2;                                     
    Noise_x=x_predic(1)-side+2*rand(1,number_returns)*side;            
    Noise_y=x_predic(3)-side+2*rand(1,number_returns)*side; 
    Noise_vx= x_predic(2)-side+rand(1,number_returns)*side;             
    Noise_vy= x_predic(4)-side+rand(1,number_returns)*side; 
    Noise=[Noise_x ;Noise_vx;Noise_y;Noise_vy]; 
     
    b=zeros(1,4); 
    b(1)=Zk(1,t); 
    b(2)=Zk(2,t); 
    b(3)=Zk(3,t); 
    b(4)=Zk(4,t); 
    y1=[Noise b'];   
    y=[]; 
    d=[]; 
    m=0; 
    
    for j=1:(number_returns+1)     
        d=y1(:,j)-Z_predic; 
        D=d'*inv(S)*d;  
        if D<=g_sigma 
            plot(y1(1,j),y1(3,j),'.');
            hold on;
            y=[y y1(:,j)];   
            m=m+1;          
        end 
    end 
     
    %% PDA_Function
    [x_putput, P] = PDA_Function(x_predic, P_predic, S, Z_predic, K,y, m);
    
    x_filter(:,t) = x_putput;
    sumx = sumx + x_filter(1,t);
    x_filter1(:,t) = x_filter(:,t);
    x_filter1(1,t) = sumx/t;
    
end

plot(target_position(1,:),target_position(3,:),'-bo');
hold on; 
plot(x_filter(1,:),x_filter(3,:),'-r*');
hold on; 
plot(x_filter1(1,:),x_filter1(3,:),'-ro');
hold off; 
axis([190 210 400 1400]);
toc

