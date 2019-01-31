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
  
simTime=300;                %仿真迭代次数
T=0.5;                     %采样时间
w2=5*2*pi/360;     %模型2转弯率3度
w3=-5*2*pi/360;    %模型3转弯率-3度
H=[1,0,0,0;0,0,1,0];                      %模型量测矩阵
G=[T^2/2,0;T,0;0,T^2/2;0,T];              %模型过程噪声加权矩阵
Q=[2^2,0;0,2^2];                                  %模型过程噪声协方差矩阵
r=5^2;                                 %20 2000
R=[r,0;0,r];                            %模型量测噪声协方差矩阵
noise_total = 6; 

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
x1=[150,20,200,-30]';  % 初始状态

% 产生量测数据
x = zeros(4,simTime);
z = zeros(2,simTime);         

x(:,1)=x0;
noise_1=H*x(:,1)+sqrt(R)*randn(2,noise_total);

x2(:,1)=x1;
noise_2=H*x2(:,1)+sqrt(R)*randn(2,noise_total);
z_noise(:,:,1) =[noise_1 noise_2];


for a=2:simTime
    if (a>=20)&&(a<=50) || (a>=120)&&(a<=180) 
        x(:,a)=F2*x(:,a-1);      
    elseif (a>=50)&&(a<=80) || (a>=210)&&(a<=270) 
        x(:,a)=F3*x(:,a-1);        
    else
        x(:,a)=F1*x(:,a-1);     
    end
    
     if (a>=10)&&(a<=60) || (a>=190)&&(a<=240) 
       x2(:,a)=F3*x2(:,a-1);      
    elseif (a>=90)&&(a<=150) || (a>=260)&&(a<=280) 
        x2(:,a)=F2*x2(:,a-1);        
    else
        x2(:,a)=F1*x2(:,a-1);     
     end
     
    noise_1=H*x(:,a)+sqrt(R)*randn(2,noise_total);
    noise_2=H*x2(:,a)+sqrt(R)*randn(2,noise_total);
    z_noise(:,:,a) =[noise_1 noise_2];
    
end

%作图部分
figure
plot(x(1,:),x(3,:),'s-'); hold on;grid on;
plot(x2(1,:),x2(3,:),'s-'); hold on;grid on;
for a=1:simTime
    plot(z_noise(1,:,a), z_noise(2,:,a),'c.');
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

%子模型参数初始化
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

%main函数
for t=1:simTime-1
    %IMM
    for index = 1:2
        x1_IMM_PDA = x1_IMM_PDA_rec(:,index) ;
        x2_IMM_PDA = x2_IMM_PDA_rec(:,index) ;
        x3_IMM_PDA = x3_IMM_PDA_rec(:,index) ;
        P1_IMM_PDA = P1_IMM_PDA_rec(:,:,index);
        P2_IMM_PDA = P2_IMM_PDA_rec(:,:,index);
        P3_IMM_PDA = P3_IMM_PDA_rec(:,:,index);
        u_IMM_PDA(:,t) = u_IMM_PDA_rec(:,t,index);   
        
        c_j(:,index)=pij'*u_IMM_PDA(:,t);
        %计算每个模型的条件概率
        ui1=(1/c_j(1,index))*pij(:,1).*u_IMM_PDA(:,t);
        ui2=(1/c_j(2,index))*pij(:,2).*u_IMM_PDA(:,t);
        ui3=(1/c_j(3,index))*pij(:,3).*u_IMM_PDA(:,t);    

        % 计算各模型滤波初始化条件
        x01(:,index)=x1_IMM_PDA*ui1(1)+x2_IMM_PDA*ui1(2)+x3_IMM_PDA*ui1(3);
        x02(:,index)=x1_IMM_PDA*ui2(1)+x2_IMM_PDA*ui2(2)+x3_IMM_PDA*ui2(3);
        x03(:,index)=x1_IMM_PDA*ui3(1)+x2_IMM_PDA*ui3(2)+x3_IMM_PDA*ui3(3);   

        %各模型滤波初始状态协方差矩阵
        P01(:,:,index)=(P1_IMM_PDA+(x1_IMM_PDA-x01(:,index))*(x1_IMM_PDA-x01(:,index))')*ui1(1)+...
                (P2_IMM_PDA+(x2_IMM_PDA-x01(:,index))*(x2_IMM_PDA-x01(:,index))')*ui1(2)+...
                (P3_IMM_PDA+(x3_IMM_PDA-x01(:,index))*(x3_IMM_PDA-x01(:,index))')*ui1(3);
        P02(:,:,index)=(P1_IMM_PDA+(x1_IMM_PDA-x02(:,index))*(x1_IMM_PDA-x02(:,index))')*ui2(1)+...
                (P2_IMM_PDA+(x2_IMM_PDA-x02(:,index))*(x2_IMM_PDA-x02(:,index))')*ui2(2)+...
                (P3_IMM_PDA+(x3_IMM_PDA-x02(:,index))*(x3_IMM_PDA-x02(:,index))')*ui2(3);
        P03(:,:,index)=(P1_IMM_PDA+(x1_IMM_PDA-x03(:,index))*(x1_IMM_PDA-x03(:,index))')*ui3(1)+...
                (P2_IMM_PDA+(x2_IMM_PDA-x03(:,index))*(x2_IMM_PDA-x03(:,index))')*ui3(2)+...
                (P3_IMM_PDA+(x3_IMM_PDA-x03(:,index))*(x3_IMM_PDA-x03(:,index))')*ui3(3); 
                   
        x_predic = F1*x01(:,index);
        P_predic = F1*P01(:,:,index)*F1'+G*Q*G';
        Z_predic(:,(index-1)*3+1) = H*x_predic;
        S(:,:,(index-1)*3+1) = H*P_predic*H'+ R;
        
        x_predic = F2*x02(:,index);
        P_predic = F2*P02(:,:,index)*F2'+G*Q*G';
        Z_predic(:,(index-1)*3+2) = H*x_predic;
        S(:,:,(index-1)*3+2) = H*P_predic*H'+ R;
        
        x_predic = F3*x03(:,index);
        P_predic = F3*P03(:,:,index)*F3'+G*Q*G';
        Z_predic(:,(index-1)*3+3) = H*x_predic;
        S(:,:,(index-1)*3+3) = H*P_predic*H'+ R;
    end
    
    %计算概率（cheap_JPDA算法）
    [y, confirm_num, U] = Cheap_JPDA_Func(6, S, Z_predic, z_noise(:,:,t+1), noise_total*2);
    
    for index = 1:2
        [x1_IMM_PDA,P1_IMM_PDA,L1] = IMMJPDA_Function(x01(:,index), P01(:,:,index), F1, H, y, confirm_num, U(:,(index-1)*3+1));
        [x2_IMM_PDA,P2_IMM_PDA,L2] = IMMJPDA_Function(x02(:,index), P02(:,:,index), F2, H, y, confirm_num, U(:,(index-1)*3+2));
        [x3_IMM_PDA,P3_IMM_PDA,L3] = IMMJPDA_Function(x03(:,index), P03(:,:,index), F3, H, y, confirm_num, U(:,(index-1)*3+3));
        %第三步--模型概率更新
        [u_IMM_PDA(:,t+1)] = Model_P_up_PDA(L1,L2,L3,c_j(:,index));
        %第四步--模型综合
        [x_pro_IMM_PDA(:,t+1),P_IMM(:,:,t+1)]=...
                                            Model_mix(x1_IMM_PDA,x2_IMM_PDA,x3_IMM_PDA,...
                                            P1_IMM_PDA,P2_IMM_PDA,P3_IMM_PDA,u_IMM_PDA(:,t));   

        x1_IMM_PDA_rec(:,index) = x1_IMM_PDA;
        x2_IMM_PDA_rec(:,index) = x2_IMM_PDA;
        x3_IMM_PDA_rec(:,index) = x3_IMM_PDA;
        P1_IMM_PDA_rec(:,:,index) = P1_IMM_PDA;
        P2_IMM_PDA_rec(:,:,index) = P2_IMM_PDA;
        P3_IMM_PDA_rec(:,:,index) = P3_IMM_PDA;
        u_IMM_PDA_rec(:,t+1,index) = u_IMM_PDA(:,t+1);    

        x_pro_IMM_PDA_rec(:,t+1,index) = x_pro_IMM_PDA(:,t+1);
    end
end

% 显示预测航迹
plot(x_pro_IMM_PDA_rec(1,:,1),x_pro_IMM_PDA_rec(3,:,1),'s-');
plot(x_pro_IMM_PDA_rec(1,:,2),x_pro_IMM_PDA_rec(3,:,2),'s-');

% % 模型概率
t=1:simTime;
figure
subplot(211)
plot(t,u_IMM_PDA_rec(1,t,1),'k.-',t,u_IMM_PDA_rec(2,t,1),'r.-',t,u_IMM_PDA_rec(3,t,1),'b.-');grid on
title('IMMJPDA Trace 1 模型概率');
axis([0 simTime-1 0 1])
xlabel('t/s'); ylabel('模型概率');
legend('模型1','模型2','模型3');
subplot(212)
plot(t,u_IMM_PDA_rec(1,t,2),'k.-',t,u_IMM_PDA_rec(2,t,2),'r.-',t,u_IMM_PDA_rec(3,t,2),'b.-');grid on
title('IMMJPDA Trace 2 模型概率');
axis([0 simTime-1 0 1])
xlabel('t/s'); ylabel('模型概率');
legend('模型1','模型2','模型3');

% 位置误差
figure
subplot(2,1,1);
t=1:simTime;
plot(t,abs(x_pro_IMM_PDA_rec(1,t,1)-x(1,t)),'LineWidth',1);hold on;grid on
plot(t,abs(x_pro_IMM_PDA_rec(1,t,2)-x2(1,t)),'LineWidth',1);hold on;grid on
title('x坐标位置跟踪误差');
xlabel('t/s'); ylabel('x-error/m');
legend('Trace 1','Trace 2');

subplot(2,1,2);
t=1:simTime;
plot(t,abs(x_pro_IMM_PDA_rec(3,t,1)-x(3,t)),'LineWidth',1);hold on;grid on
plot(t,abs(x_pro_IMM_PDA_rec(3,t,2)-x2(3,t)),'LineWidth',1);hold on;grid on
title('y坐标位置跟踪误差');
xlabel('t/s'); ylabel('y-error/m');
legend('Trace 1','Trace 2');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% 相关子函数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%二、模型概率更新函数
function [u]=Model_P_up(r1,r2,r3,S1,S2,S3,c_j)
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

function [y,confirm_num,U] = Cheap_JPDA_Func(trace_num, S, Z_predic, z_noise, measure_num)
    y = [];
    confirm_num = 0;
    %生成关联矩阵
    for j=1:measure_num 
        flag=0;
        for i=1:trace_num
            d=z_noise(:,j)-Z_predic(:,i);
            D=d'*inv(S(:,:,i))*d; 
            if D<=10                                                    
               flag=1;
            end
        end
        if flag==1   
           y=[y z_noise(:,j)];                                                     
           confirm_num=confirm_num+1;                                                
        end
    end 
        
    Gjt = zeros(trace_num,confirm_num);
    for i = 1:trace_num
        temp = 0;
        for j = 1:confirm_num
            v1 = y(:,j)-Z_predic(:,i);
            v = (-0.5*v1'*inv(S(:,:,i))*v1);
            Gjt(j,i) = exp(v)/(sqrt(2*3.14*det(S(:,:,i))));%似然函数
            temp = temp + Gjt(j,i);
        end
        St(i) =  temp;
    end
    
    B = 0;
    for i = 1:trace_num
        ST = St(i);
        for j = 1:confirm_num
            GJT = Gjt(j,i);
            Sj = sum(Gjt(j,1:trace_num));        
            U(j,i) = GJT/((ST+Sj-GJT+B));
        end
    end
      
    for i = 1:trace_num
        U(j+1,i) = 1 - sum(U(1:j,i));
    end
end

