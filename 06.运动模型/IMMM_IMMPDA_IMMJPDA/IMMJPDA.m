% 基于IMM算法的目标跟踪
%   20181228: 1）更新IMMJPDA算法
clc;clear all;close all;

tic
n=50;
global Pd Pg gamma G Q R noise_total;
Pd=1;       
Pg=0.97;   
gamma = 1; 

g_sigma = 1;     
  
simTime=100;      %仿真迭代次数
T=0.5;                     %采样时间
w2=5*2*pi/360;     %模型2转弯率3度
w3=-5*2*pi/360;    %模型3转弯率-3度
H=[1,0,0,0;0,0,1,0];                      %模型量测矩阵
G=[T^2/2,0;T,0;0,T^2/2;0,T];              %模型过程噪声加权矩阵
Q=[2^2,0;0,2^2];                                  %模型过程噪声协方差矩阵
r=3^2;                                 %20 2000
R=[r,0;0,r];                            %模型量测噪声协方差矩阵

F1=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];     %模型1状态转移矩阵

F2=[1,sin(w2*T)/w2,0,(cos(w2*T)-1)/w2;
    0,cos(w2*T),0,sin(w2*T);
    0,(1-cos(w2*T))/w2,1,sin(w2*T)/w2;
    0,-sin(w2*T),0,cos(w2*T)];            %模型2状态转移矩阵 左转弯

F3=[1,sin(w3*T)/w3,0,(cos(w3*T)-1)/w3;
    0,cos(w3*T),0,sin(w3*T);
    0,(1-cos(w3*T))/w3,1,sin(w3*T)/w3;
    0,-sin(w3*T),0,cos(w3*T)];            %模型3状态转移矩阵 右转弯

x0=[100,20,100,20]';  % 初始状态
x1=[400,20,600,-30]';  % 初始状态

% 产生量测数据
x = zeros(4,simTime);
z = zeros(2,simTime);         

x(:,1)=x0;
z(:,1)=H*x(:,1)+sqrt(R)*randn(2,1);

x2(:,1)=x1;
z2(:,1)=H*x2(:,1)+sqrt(R)*randn(2,1);


noise_total = 6; 
z_noise = zeros((noise_total+1)*2,2,simTime);   
ellipse_Volume= pi*gamma*30;
side=sqrt((ellipse_Volume*gamma+1)/gamma)/2;
Noise_x= z(1,1)+side-2*rand(1,noise_total)*side; 
Noise_y= z(2,1)+side-2*rand(1,noise_total)*side;
z_noise_1 =[[Noise_x;Noise_y]'; z(:,1)'];
Noise_x= z2(1,1)+side-2*rand(1,noise_total)*side; 
Noise_y= z2(2,1)+side-2*rand(1,noise_total)*side;  
z_noise(:,:,1) =[z_noise_1;[Noise_x;Noise_y]'; z2(:,1)'];

for a=2:simTime
    if (a>=20)&&(a<=50) 
        x(:,a)=F2*x(:,a-1);      
    elseif (a>=50)&&(a<=80) 
        x(:,a)=F3*x(:,a-1);        
    else
        x(:,a)=F1*x(:,a-1);     
    end
    
     if (a>=1)&&(a<=30) 
       x2(:,a)=F2*x2(:,a-1);      
    elseif (a>=30)&&(a<=90) 
        x2(:,a)=F3*x2(:,a-1);        
    else
        x2(:,a)=F1*x2(:,a-1);     
    end
    
    z(:,a)=H*x(:,a)+sqrt(R)*randn(2,1);
    z2(:,a)=H*x2(:,a)+sqrt(R)*randn(2,1);
    
    Noise_x= z(1,a)+side-2*rand(1,noise_total)*side; 
    Noise_y= z(2,a)+side-2*rand(1,noise_total)*side;
    z_noise_1 =[[Noise_x;Noise_y]'; z(:,a)'];
    Noise_x= z2(1,a)+side-2*rand(1,noise_total)*side; 
    Noise_y= z2(2,a)+side-2*rand(1,noise_total)*side;  
    z_noise(:,:,a) =[z_noise_1;[Noise_x;Noise_y]'; z2(:,a)'];
end

%作图部分
figure
plot(x(1,:),x(3,:),'b*-'); hold on;grid on;
plot(x2(1,:),x2(3,:),'r*-'); hold on;grid on;
for a=1:simTime
    plot(z_noise(:,1,a), z_noise(:,2,a),'k.');
end

x_pro_IMM(:,1)=x0;
x_pro_IMM_PDA(:,1)=x0;
%模型转移概率矩阵
pij=[0.95,0.025,0.025;
       0.025,0.95,0.025;
       0.025,0.025,0.95];    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%基本IMM算法参数
u_IMM=zeros(3,simTime);
%IMM算法模型概率
u_IMM(:,1)=[0.8,0.1,0.1]';  
%IMM算法各模型初始状态
x1_IMM=x0;x2_IMM=x0;x3_IMM=x0; 

 %初始状态协方差矩阵
P0=diag([10,10,10,10]); 
P1_IMM=P0;P2_IMM=P0;P3_IMM=P0;
P_IMM(:,:,1)=P0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMMPDA
u_IMM_PDA=zeros(3,simTime);
u_IMM_PDA(:,1)=[0.8,0.1,0.1]';  
x1_IMM_PDA=x0;x2_IMM_PDA=x0;x3_IMM_PDA=x0; 

P0=diag([10,10,10,10]);  
P1_IMM_PDA=P0;P2_IMM_PDA=P0;P3_IMM_PDA=P0;
P_IMM_PDA(:,:,1)=P0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
x4 = x1_IMM;
x4_rec(:,1) = x4;
P4 = P0;

x1_IMM_PDA_rec(:,1) = x0;
x2_IMM_PDA_rec(:,1) = x0;
x3_IMM_PDA_rec(:,1) = x0;
P1_IMM_PDA_rec(:,:,1) = P0;
P2_IMM_PDA_rec(:,:,1) = P0;
P3_IMM_PDA_rec(:,:,1) = P0;
u_IMM_PDA_rec(:,1,1) = u_IMM_PDA(:,1);    

x1_IMM_PDA_rec(:,2) = x1;
x2_IMM_PDA_rec(:,2) = x1;
x3_IMM_PDA_rec(:,2) = x1;
P1_IMM_PDA_rec(:,:,2) = P0;
P2_IMM_PDA_rec(:,:,2) = P0;
P3_IMM_PDA_rec(:,:,2) = P0;
u_IMM_PDA_rec(:,1,2) = u_IMM_PDA(:,1);    

x_pro_IMM_PDA_rec = zeros(4,100,2);
x_pro_IMM_PDA_rec(:,1,1) = x0;
x_pro_IMM_PDA_rec(:,1,2) = x1;

for t=1:simTime-1
  
    %IMMPDA算法
    for jjj = 1:2
        x1_IMM_PDA = x1_IMM_PDA_rec(:,jjj) ;
        x2_IMM_PDA = x2_IMM_PDA_rec(:,jjj) ;
        x3_IMM_PDA = x3_IMM_PDA_rec(:,jjj) ;
        P1_IMM_PDA = P1_IMM_PDA_rec(:,:,jjj);
        P2_IMM_PDA = P2_IMM_PDA_rec(:,:,jjj);
        P3_IMM_PDA = P3_IMM_PDA_rec(:,:,jjj);
        u_IMM_PDA(:,t) = u_IMM_PDA_rec(:,t,jjj);   
        
        c_j=pij'*u_IMM_PDA(:,t);
        %计算每个模型的条件概率
        ui1=(1/c_j(1))*pij(:,1).*u_IMM_PDA(:,t);
        ui2=(1/c_j(2))*pij(:,2).*u_IMM_PDA(:,t);
        ui3=(1/c_j(3))*pij(:,3).*u_IMM_PDA(:,t);    

        % 计算各模型滤波初始化条件
        x01=x1_IMM_PDA*ui1(1)+x2_IMM_PDA*ui1(2)+x3_IMM_PDA*ui1(3);
        x02=x1_IMM_PDA*ui2(1)+x2_IMM_PDA*ui2(2)+x3_IMM_PDA*ui2(3);
        x03=x1_IMM_PDA*ui3(1)+x2_IMM_PDA*ui3(2)+x3_IMM_PDA*ui3(3);   

        %各模型滤波初始状态协方差矩阵
        P01=(P1_IMM_PDA+(x1_IMM_PDA-x01)*(x1_IMM_PDA-x01)')*ui1(1)+...
                (P2_IMM_PDA+(x2_IMM_PDA-x01)*(x2_IMM_PDA-x01)')*ui1(2)+...
                (P3_IMM_PDA+(x3_IMM_PDA-x01)*(x3_IMM_PDA-x01)')*ui1(3);
        P02=(P1_IMM_PDA+(x1_IMM_PDA-x02)*(x1_IMM_PDA-x02)')*ui2(1)+...
                (P2_IMM_PDA+(x2_IMM_PDA-x02)*(x2_IMM_PDA-x02)')*ui2(2)+...
                (P3_IMM_PDA+(x3_IMM_PDA-x02)*(x3_IMM_PDA-x02)')*ui2(3);
        P03=(P1_IMM_PDA+(x1_IMM_PDA-x03)*(x1_IMM_PDA-x03)')*ui3(1)+...
                (P2_IMM_PDA+(x2_IMM_PDA-x03)*(x2_IMM_PDA-x03)')*ui3(2)+...
                (P3_IMM_PDA+(x3_IMM_PDA-x03)*(x3_IMM_PDA-x03)')*ui3(3); 
                   
        x_predic = F1*x01;
        P_predic = F1*P01*F1'+G*Q*G';
        Z_predic(:,(jjj-1)*3+1) = H*x_predic;
        S(:,:,(jjj-1)*3+1) = H*P_predic*H'+ R;
        
        x_predic = F2*x02;
        P_predic = F2*P02*F2'+G*Q*G';
        Z_predic(:,(jjj-1)*3+2) = H*x_predic;
        S(:,:,(jjj-1)*3+2) = H*P_predic*H'+ R;
        
        x_predic = F3*x03;
        P_predic = F3*P03*F3'+G*Q*G';
        Z_predic(:,(jjj-1)*3+3) = H*x_predic;
        S(:,:,(jjj-1)*3+3) = H*P_predic*H'+ R;
    end
    
    c = 6;
    n2 = noise_total*2+2;
    y = [];
    m1 = 0;
    %生成关联矩阵
    for j=1:n2 
        flag=0;
        for i=1:c
            d=z_noise(j,:,t+1)'-Z_predic(:,i);
            D=d'*inv(S(:,:,i))*d; 
            if D<=10                                                    
               flag=1;
            end
        end
        if flag==1   
           y=[y z_noise(j,:,t+1)'];                                                     
           m1=m1+1;                                                
        end
    end 
        
    Gjt = zeros(c,m1);
    for i = 1:c
        temp = 0;
        for j = 1:m1
            v1 = y(:,j)-Z_predic(:,i);
            v = (-0.5*v1'*inv(S(:,:,i))*v1);
            Gjt(j,i) = exp(v)/(sqrt(2*3.14*det(S(:,:,i))));%似然函数
            temp = temp + Gjt(j,i);
        end
        St(i) =  temp;
    end
    
    B = 0;
    for i = 1:c
        ST = St(i);
        for j = 1:m1
            GJT = Gjt(j,i);
            Sj = sum(Gjt(j,1:c));        
            U(j,i) = GJT/((ST+Sj-GJT+B));
        end
    end
      
    for i = 1:c
        U(j+1,i) = 1 - sum(U(1:j,i));
    end
    
    for jjj = 1:2
        x1_IMM_PDA = x1_IMM_PDA_rec(:,jjj) ;
        x2_IMM_PDA = x2_IMM_PDA_rec(:,jjj) ;
        x3_IMM_PDA = x3_IMM_PDA_rec(:,jjj) ;
        P1_IMM_PDA = P1_IMM_PDA_rec(:,:,jjj)  ;
        P2_IMM_PDA = P2_IMM_PDA_rec(:,:,jjj)  ;
        P3_IMM_PDA = P3_IMM_PDA_rec(:,:,jjj)  ;
        u_IMM_PDA(:,t) = u_IMM_PDA_rec(:,t,jjj);   
        
        c_j=pij'*u_IMM_PDA(:,t);
        %计算每个模型的条件概率
        ui1=(1/c_j(1))*pij(:,1).*u_IMM_PDA(:,t);
        ui2=(1/c_j(2))*pij(:,2).*u_IMM_PDA(:,t);
        ui3=(1/c_j(3))*pij(:,3).*u_IMM_PDA(:,t);    

        % 计算各模型滤波初始化条件
        x01=x1_IMM_PDA*ui1(1)+x2_IMM_PDA*ui1(2)+x3_IMM_PDA*ui1(3);
        x02=x1_IMM_PDA*ui2(1)+x2_IMM_PDA*ui2(2)+x3_IMM_PDA*ui2(3);
        x03=x1_IMM_PDA*ui3(1)+x2_IMM_PDA*ui3(2)+x3_IMM_PDA*ui3(3);   

        %各模型滤波初始状态协方差矩阵
        P01=(P1_IMM_PDA+(x1_IMM_PDA-x01)*(x1_IMM_PDA-x01)')*ui1(1)+...
                (P2_IMM_PDA+(x2_IMM_PDA-x01)*(x2_IMM_PDA-x01)')*ui1(2)+...
                (P3_IMM_PDA+(x3_IMM_PDA-x01)*(x3_IMM_PDA-x01)')*ui1(3);
        P02=(P1_IMM_PDA+(x1_IMM_PDA-x02)*(x1_IMM_PDA-x02)')*ui2(1)+...
                (P2_IMM_PDA+(x2_IMM_PDA-x02)*(x2_IMM_PDA-x02)')*ui2(2)+...
                (P3_IMM_PDA+(x3_IMM_PDA-x02)*(x3_IMM_PDA-x02)')*ui2(3);
        P03=(P1_IMM_PDA+(x1_IMM_PDA-x03)*(x1_IMM_PDA-x03)')*ui3(1)+...
                (P2_IMM_PDA+(x2_IMM_PDA-x03)*(x2_IMM_PDA-x03)')*ui3(2)+...
                (P3_IMM_PDA+(x3_IMM_PDA-x03)*(x3_IMM_PDA-x03)')*ui3(3); 
            
        [x1_IMM_PDA,P1_IMM_PDA,L1] = IMMJPDA_Function(x01, P01, F1, H, y, m1, U(:,(jjj-1)*3+1));
        [x2_IMM_PDA,P2_IMM_PDA,L2] = IMMJPDA_Function(x02, P02, F2, H, y, m1, U(:,(jjj-1)*3+2));
        [x3_IMM_PDA,P3_IMM_PDA,L3] = IMMJPDA_Function(x03, P03, F3, H, y, m1, U(:,(jjj-1)*3+3));
        %第三步--模型概率更新
        [u_IMM_PDA(:,t+1)] = Model_P_up_PDA(L1,L2,L3,c_j);
        %第四步--模型综合
        [x_pro_IMM_PDA(:,t+1),P_IMM(:,:,t+1)]=...
                                            Model_mix(x1_IMM_PDA,x2_IMM_PDA,x3_IMM_PDA,...
                                            P1_IMM_PDA,P2_IMM_PDA,P3_IMM_PDA,u_IMM_PDA(:,t));   

        x1_IMM_PDA_rec(:,jjj) = x1_IMM_PDA;
        x2_IMM_PDA_rec(:,jjj) = x2_IMM_PDA;
        x3_IMM_PDA_rec(:,jjj) = x3_IMM_PDA;
        P1_IMM_PDA_rec(:,:,jjj) = P1_IMM_PDA;
        P2_IMM_PDA_rec(:,:,jjj) = P2_IMM_PDA;
        P3_IMM_PDA_rec(:,:,jjj) = P3_IMM_PDA;
        u_IMM_PDA_rec(:,t+1,jjj) = u_IMM_PDA(:,t+1);    

        x_pro_IMM_PDA_rec(:,t+1,jjj) = x_pro_IMM_PDA(:,t+1);
    end
end


% figure
% plot(x(1,:),x(3,:),'r*-'); hold on;grid on;
% plot(x2(1,:),x2(3,:),'b*-'); hold on;grid on;
plot(x_pro_IMM_PDA_rec(1,:,1),x_pro_IMM_PDA_rec(3,:,1),'rs-');
plot(x_pro_IMM_PDA_rec(1,:,2),x_pro_IMM_PDA_rec(3,:,2),'bs-');
% plot(x4_rec(1,:),x4_rec(3,:),'s-');
% legend('Real','IMM','IMMPDA','CT');

% % 模型概率
t=1:simTime;
figure
subplot(121)
plot(t,u_IMM_PDA_rec(1,t,1),'k.-',t,u_IMM_PDA_rec(2,t,1),'r.-',t,u_IMM_PDA_rec(3,t,1),'b.-');grid on
title('IMM算法模型概率曲线');
axis([0 simTime-1 0 1])
xlabel('t/s'); ylabel('模型概率');
legend('模型1','模型2','模型3');
subplot(122)
plot(t,u_IMM_PDA_rec(1,t,2),'k.-',t,u_IMM_PDA_rec(2,t,2),'r.-',t,u_IMM_PDA_rec(3,t,2),'b.-');grid on
title('IMM算法模型概率曲线');
axis([0 simTime-1 0 1])
xlabel('t/s'); ylabel('模型概率');
legend('模型1','模型2','模型3');

% 位置误差
% figure
% subplot(2,1,1);
% t=1:simTime;
% plot(t,abs(x_pro_IMM(1,t)-x(1,t)),'LineWidth',1);hold on;grid on
% plot(t,abs(x_pro_IMM_PDA(1,t)-x(1,t)),'LineWidth',1);hold on;grid on
% plot(t,abs(x4_rec(1,t)-x(1,t)),'LineWidth',1);grid on;
% title('x坐标位置跟踪误差');
% xlabel('t/s'); ylabel('x-error/m');
% legend('IMM','IMMPDA','CT');
% 
% subplot(2,1,2);
% t=1:simTime;
% plot(t,abs(x_pro_IMM(3,t)-x(3,t)),'LineWidth',1);hold on;grid on
% plot(t,abs(x_pro_IMM_PDA(3,t)-x(3,t)),'LineWidth',1);hold on;grid on
% plot(t,abs(x4_rec(3,t)-x(3,t)),'LineWidth',1);grid on;
% title('y坐标位置跟踪误差');
% xlabel('t/s'); ylabel('y-error/m');
% legend('IMM','IMMPDA','CT');


%二、模型概率更新函数
function [u]=Model_P_up(r1,r2,r3,S1,S2,S3,c_j)
%模型概率更新函数
%计算似然函数
Lfun1=(1/sqrt(abs(2*pi*(det(S1)))))*exp((-1/2)*(r1'*inv(S1)*r1));   %Lfun1=1/(r1'*inv(S1)*r1);
Lfun2=(1/sqrt(abs(2*pi*(det(S2)))))*exp((-1/2)*(r2'*inv(S2)*r2));   %Lfun2=1/(r2'*inv(S2)*r2);
Lfun3=(1/sqrt(abs(2*pi*(det(S3)))))*exp((-1/2)*(r3'*inv(S3)*r3));    %Lfun3=1/(r3'*inv(S3)*r3);
c=[Lfun1,Lfun2,Lfun3]'.*c_j;
% 再归一化
u=[Lfun1,Lfun2,Lfun3]'.*c_j*(1/sum(c));

end

%二、模型概率更新函数
function [u]=Model_P_up_PDA(Lfun1,Lfun2,Lfun3,c_j)
%模型概率更新函数
c=[Lfun1,Lfun2,Lfun3]'.*c_j;
% 再归一化
u=[Lfun1,Lfun2,Lfun3]'.*c_j*(1/sum(c));

end

%三、状态混合函数
function [x_pro,P]=Model_mix(x1,x2,x3,P1,P2,P3,u)
x_pro=x1*u(1)+x2*u(2)+x3*u(3);
P=(P1+(x1-x_pro)*(x1-x_pro)')*u(1)+...
      (P2+(x2-x_pro)*(x2-x_pro)')*u(2)+...
      (P3+(x3-x_pro)*(x3-x_pro)')*u(3);
end

%四、Kalman滤波函数
function [X,P,e,S]=Kalman(X_Forward,P_Forward,z_noise,F,G,Q,H,R)
%卡尔曼滤波2012.2.27   IMM专用，参数略有不同
    x_predic = F*X_Forward;
    P_predic = F*P_Forward*F'+G*Q*G'; 
    Z_predic = H*x_predic;
    S = H*P_predic*H'+ R; 
    K = P_predic*H'*inv(S);
    y=[];m=0;
    for j=1:size(z_noise,1)
        d=z_noise(j,:)'-Z_predic;  
        D(j)=d'*inv(S)*d;  
    end
    [~,index] = min(D);
    Z = z_noise(index,:)';
    e = Z-Z_predic; %残差

    % 修正滤波值和误差协方差阵
    X=x_predic+K*e;
    M=K*H;
    n=size(M);
    I=eye(n);
    P=(I-K*H)*P_predic*(I-K*H)'+ K*R*K';
end

function [X_out, P_out] = cheap_JPDA_xxx(x_predic_in, P_in ,Z_predic_in ,S_in ,num, y, m, ffun,hfun,Q)
    Gjt = zeros(m,num);
    for i = 1:num
        temp = 0;
        for j = 1:m
            v1 = y(:,j)-Z_predic_in(:,i);
            v = (-0.5*v1'*inv(S_in(:,:,i))*v1);
            Gjt(j,i) = exp(v)/(2*3.14*sqrt(det(S_in(:,:,i))));
            temp = temp + Gjt(j,i);
        end
        St(i) =  temp;
    end
    
    B = 0;
    for i = 1:num
        ST = St(i);
        for j = 1:m
            GJT = Gjt(j,i);
            Sj = 0;
            for k =1:num
                Sj = Sj + Gjt(j,k);
            end        
            U(j,i) = GJT/((ST+Sj-GJT+B));
        end
    end
      
    for i = 1:num
        U(j+1,i) = 1 - sum(U(1:j,i));
    end
    
    %Kalman filter
    for i=1:num                                                                
        P_predic = ffun*P_in(:,:,i)*ffun'+Q;
        K(:,:,i)= P_predic*hfun'*inv(S_in(:,:,i));
        P(:,:,i)= P_predic-(1-U(m+1,i))*K(:,:,i)*S_in(:,:,i)*K(:,:,i)';
    end

    for i=1:num
        a=0;         
        b=0;

        x_filter_temp=0;
        for j=1:m 
            x_filter_temp=x_filter_temp+U(j,i)*(x_predic_in(:,i)+ K (:,:,i)*(y(:,j)- Z_predic_in(:,i)));
        end
        x_filter_temp=U(j+1,i)*x_predic_in(:,i)+x_filter_temp;
        x_filter(:,i)=x_filter_temp;

        for j=1:m+1
            if j==m+1
                a=x_predic_in(:,i);
            else
               a=x_predic_in(:,i)+ K (:,:,i)*(y(:,j)- Z_predic_in(:,i));
            end
            b=b+U(j,i)*(a*a'-x_filter_temp*x_filter_temp');
        end
        P_out(:,:,i)=P(:,:,i)+b; 
        X_out(:,i)=x_filter(:,i);
    end
       
end


