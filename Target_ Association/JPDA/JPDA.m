clc;
clear all;
close all;

global Pd Pg g_sigma lambda c gamma;
Pd=1;         %检测概率
Pg=0.99;      %正确量测落入跟踪门内得概率
g_sigma=9.21; %门限（卡方值）
lambda=8;     %杂波密度
c=2;          %目标个数
gamma=lambda*10^(-6); 

n=50;   %采样次数
T=1;    %T为采样间隔
MC_number=5;     %Monte Carlo仿真次数
                                                                             
target_position=[1500 300 800 400; 1000 400 3000 300];        %目标的起始位置和速度(m,m/s)                   
JPDAF(target_position,n,T,MC_number,c);     


function  JPDAF(target_position,n,T,MC_number,c)
global g_sigma gamma;                                                             %每一个单位面积(km^2)内产生lambda个杂波
% Target_measurement=zeros(c,2,n);                                                   %目标观测互联存储矩阵
target_delta=[100 100 100];                                                            %目标对应的观测标准差                    
P=zeros(4,4,c);                                                                    %协方差矩阵
P1=[target_delta(1)^2 0 0 0;0 0.01 0 0;0 0 target_delta(1)^2 0;0 0 0 0.01];        %初始协方差矩阵 
P(:,:,1)=P1;
P(:,:,2)=P1;
A = [1 T 0 0;
    0 1 0 0;
    0 0 1 T;
    0 0 0 1];   %状态转移矩阵
C = [1 0 0 0;
    0 0 1 0];  %观测矩阵

R=[target_delta(1)^2 0;0 target_delta(1)^2];                                       %观测协方差矩阵
Q=[4 0;0 4];                                                                       %系统过程噪声协方差
G=[T^2/2 0;T 0;0 T^2/2;0 T];                                                       %过程噪声矩阵
x_filter=zeros(4,c,n);                                                             %存储目标的各时刻的滤波值
x_filter1=zeros(4,c,n,MC_number);                                                  %MC_number次Montle Carlo仿真所得全部结果存储矩阵
data_measurement=zeros(c,2,n);                                                     %观测存储矩阵
data_measurement1=zeros(c,4,n);                                                    %实际位置坐标x,y矩阵   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  产生目标的实际位置  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_measurement1(:,:,1)=target_position;                                          %实际位置矩阵初始化 
for i=1:c
    for ii=2:n    
        %计算两个目标的实际位置 
        data_measurement1(i,:,ii)=(A*data_measurement1(i,:,ii-1)')'+(G*sqrt(Q)*(randn(2,1)))';        
    end
end

for M=1:MC_number
%%%%%%%%%%%%%%%%%%%%    
%%%  1.产生路径  %%%
%%%%%%%%%%%%%%%%%%%%
Noise=[];
for i=1:n
    for j=1:c                                                                      
        %得到两个目标的观测位置（含噪声）
        data_measurement(j,1,i)=data_measurement1(j,1,i)+rand(1)*target_delta(j);
        data_measurement(j,2,i)=data_measurement1(j,3,i)+rand(1)*target_delta(j); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  2.产生杂波,并确定有效观测  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=zeros(2,2,c);

%存储两个目标的观测预测值,即只包括x,y坐标
Z_predic=zeros(2,2);                                                               

%存储两个目标的状态预测值,即包括x,y坐标和x,y方向速度
x_predic=zeros(4,2);                                                               

ellipse_Volume=zeros(1,2);

 %存储目标1的杂波
NOISE_sum_a=[]; 
 %存储目标2的杂波
NOISE_sum_b=[];
 %存储目标3的杂波
NOISE_sum_c=[];

for t=1:n
    y1=[];
    y=[];
    Noise=[];
    NOISE=[];
    
    for i=1:c      
        if t~=1
            %用前一时刻的滤波值来预测当前的值(kalman滤波的第一个表达式)
            x_predic(:,i) = A*x_filter(:,i,t-1);      
%             x_predic(:,i)=data_measurement1(i,:,t)'; 
        else
            %第一次采样我们用真实位置当预测值
            x_predic(:,i)=target_position(i,:)';                                                   
        end
        
        %更新x_predic协方差矩阵(kalman滤波的第二个表达式) 
        P_predic=A*P(:,:,i)*A'+G*Q*G';                                                
        
        %提取预测值的x,y坐标，舍弃x,y速度                                  
        Z_predic(:,i)=C*x_predic(:,i);                                                 
        R=[target_delta(i)^2 0; 0 target_delta(i)^2];
    
        %S为新息协方差
        S(:,:,i)=C*P_predic*C'+R;   
        
        %下列步骤是为了生成多余回波，可以理解为虚假目标
        %量测确认面积区域
        ellipse_Volume(i)= pi*g_sigma*sqrt(det(S(:,:,i)));          %计算椭圆跟踪门的面积   
%          number_returns=floor(ellipse_Volume(i)*gamma+1);           %椭圆跟踪门内的错误回波数
        number_returns = 2;       
        side=sqrt((ellipse_Volume(i)*gamma+1)/gamma)/2;            %将椭圆跟踪门等效为正方形，并求出正方形边长的二分之一
        Noise_x=x_predic(1,i)+side-2*rand(1,number_returns)*side;  %在预测值周围产生多余回波。注意：当某一次number_returns小于等于0时会出错，再运行一次即可。
        Noise_y=x_predic(3,i)+side-2*rand(1,number_returns)*side;    
        Noise=[Noise_x;Noise_y];
        NOISE=[NOISE Noise];
        
        if i==1
            NOISE_sum_a =[NOISE_sum_a Noise]; 
        end
        
       if i==2
            NOISE_sum_b =[NOISE_sum_b Noise]; 
       end
       
    end
    
    b=zeros(1,2);
    b(1)=data_measurement(1,1,t);
    b(2)=data_measurement(1,2,t);
    y1=[NOISE b'];                %将接收到的所有的回波存在y1中,包括杂波和观测
    b(1)=data_measurement(2,1,t);
    b(2)=data_measurement(2,2,t);
    y1=[y1 b'];                  
    [n1,n2]=size(y1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  3.产生观测确认矩阵Q2  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m1=0;                         %记录有效观测个数
    [n1,n2]=size(y1);
    Q1=zeros(100,3);
    for j=1:n2 
        flag=0;
        for i=1:c
            %计算量测值和预测值的残差
            d=y1(:,j)-Z_predic(:,i);
            D=d'*inv(S(:,:,i))*d; 
            %判断是否满足门限要求（这里其实是按照卡方分布来求概率）
            if D<=g_sigma                                                    
               flag=1;
               Q1(m1+1,1)=1;
               Q1(m1+1,i+1)=1;
            end
        end
        if flag==1   
           y=[y y1(:,j)];                                                      %把落入跟踪门中的所有回波放入y中
           m1=m1+1;                                                            %记录有效观测个数
        end
    end
    Q2=Q1(1:m1,1:3);
    [U] = JPDA_func(Q2,m1,ellipse_Volume,y,Z_predic,S);

%%  7.Kalman滤波开始  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:c                                                                          %更新协方差矩阵
    P_predic = A*P(:,:,i)*A'+G*Q*G';
    %卡尔曼增益
    K(:,:,i)= P_predic*C'*inv(S(:,:,i));
    %详见《雷达数据。。。》P165(8.130部分1)
    P(:,:,i)= P_predic-(1-U(m1+1,i))*K(:,:,i)*S(:,:,i)*K(:,:,i)';
end

%% 
for i=1:c
    a=0;         
    b=0;
    
    %% 计算状态预测值
    x_filter2=0;
    for j=1:m1 
        x_filter2=x_filter2+U(j,i)*(x_predic(:,i)+ K (:,:,i)*(y(:,j)- Z_predic(:,i)));
    end
    x_filter2=U(j+1,i)*x_predic(:,i)+x_filter2;
    x_filter(:,i,t)=x_filter2;
    
    %% 计算状态协方差
    for j=1:m1+1
        if j==m1+1
            a=x_predic(:,i);
        else
            %量测预测值
           a=x_predic(:,i)+ K (:,:,i)*(y(:,j)- Z_predic(:,i));
        end
        %详见《雷达数据。。。》P135 公式7.102第二部分
        b=b+U(j,i)*(a*a'-x_filter2*x_filter2');
    end
    %协方差
    %详见《雷达数据。。。》P136 公式7.102
    P(:,:,i)=P(:,:,i)+b; 
    x_filter1(:,i,t,M)=x_filter(:,i,t);
end
end
end
 

%% 做图部分
% x_filter=sum(x_filter1,4)/MC_number;   %滤波值作平均
%%%%%%%%%%%%%%%%%%%%
%%%  1.滤波结果  %%%
%%%%%%%%%%%%%%%%%%%%
figure;
%目标a,b的观测位置
for i=1:c
    a=zeros(1,n);
    b=zeros(1,n);
    for j=1:n
        a(j)=data_measurement(i,1,j);
        b(j)=data_measurement(i,2,j);
    end
    if i==1
       plot(a(:),b(:),'r*')
    end
    if i ==2 
       plot(a(:),b(:),'b*')
    end
    hold on;
end
%目标a,b的杂波位置
plot(NOISE_sum_a(1,:),NOISE_sum_a(2,:),'r.');
plot(NOISE_sum_b(1,:),NOISE_sum_b(2,:),'b.');

hold on;
%目标a,b的估计位置
a(1:n) = x_filter(1,1,1:n);
b(1:n) = x_filter(3,1,1:n);
plot(a,b,'r-');
a(1:n) = x_filter(1,2,1:n);
b(1:n) = x_filter(3,2,1:n);
plot(a,b,'b-');
legend('目标a的观测位置','目标b的观测位置','目标a的杂波','目标b的杂波','目标a的估计位置','目标b的估计位置');grid;


%%%%%%%%%%%%%%%%%%%%
%%%  2.速度误差  %%%
%%%%%%%%%%%%%%%%%%%%
figure;
a=0;
c1=zeros(c,n);
for j=1:n
    for i=1:MC_number    %最小均方误差
        a=(x_filter1(1,1,j,i)-data_measurement1(1,1,j))^2+(x_filter1(3,1,j,i)-data_measurement1(1,3,j))^2;
        c1(1,j)=c1(1,j)+a;
    end
        c1(1,j)=sqrt(c1(1,j)/MC_number);
end

temp=c1(1,:);
a_extra=zeros(2,n);
b_extra=zeros(1,n);
c_extra=zeros(1,n);
a_extra(1,:)=temp;
a_extra(2,:)=1:1:n;
b_extra=a_extra(1,:);
[c_extra,pos]=sort(b_extra);                                                       %pos为排序后的下标,c为第一行的排序结果;
a_extra(2,:)=a_extra(2,pos);                                                       %第二行按照第一行排序的下标对应
a_extra(1,:)=c_extra;                                                              %第一行结果重新赋给a 的第一行;
str1=num2str(a_extra(2,n));
str2=num2str(a_extra(1,n));
str=strcat('\itN=',str1,'\itError=',str2,'(m)');
text(a_extra(2,n),0.8*a_extra(1,n),str);
hold on;
plot([a_extra(2,n) a_extra(2,n)],[0 a_extra(1,n)],'r');
hold on;
plot(1:n,c1(1,:),'r:'); 
hold on;
a=0;

for j=1:n
    for i=1:MC_number                                                              %最小均方误差
        a=(x_filter1(1,2,j,i)-data_measurement1(2,1,j))^2+(x_filter1(3,2,j,i)-data_measurement1(2,3,j))^2;
        c1(2,j)=c1(2,j)+a;
    end
        c1(2,j)=sqrt(c1(2,j)/MC_number);
end

temp=c1(2,:);
a_extra=zeros(2,n);
b_extra=zeros(1,n);
c_extra=zeros(1,n);
a_extra(1,:)=temp;
a_extra(2,:)=1:1:n;
b_extra=a_extra(1,:);
[c_extra,pos]=sort(b_extra);                                                       %pos为排序后的下标,c为第一行的排序结果;
a_extra(2,:)=a_extra(2,pos);                                                       %第二行按照第一行排序的下标对应
a_extra(1,:)=c_extra;                                                              %第一行结果重新赋给a 的第一行;
str1=num2str(a_extra(2,n));
str2=num2str(a_extra(1,n));
str=strcat('\itN=',str1,'\itError=',str2,'(m)');
text(a_extra(2,n),0.8*a_extra(1,n),str);
hold on;
plot([a_extra(2,n) a_extra(2,n)],[0 a_extra(1,n)],'b');
hold on;
plot(1:n,c1(2,:),'b:');

xlabel('times'),ylabel('测量值与估计值均方差/m');
legend('目标a的误差最大值','目标a的误差','目标b的误差最大值','目标b的误差');grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Revised on 26th June 2008 by wangzexun  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
extra11=zeros(1,n);
extra12=zeros(1,n);
extra13=zeros(1,n);
for j=1:n
    extra11(1,j)=sqrt(x_filter(1,1,j)-data_measurement1(1,1,j))^2+(x_filter(3,1,j)-data_measurement1(1,3,j))^2;
    extra12(1,j)=sqrt((data_measurement(1,1,j)-data_measurement1(1,1,j))^2+(data_measurement(1,2,j)-data_measurement1(1,3,j))^2);
    extra13(1,j)=extra12(1,j)/extra11(1,j);
end
plot(1:n,extra13(1,:),'k:'); 
xlabel('times'),ylabel('RMSE of a');
grid;

figure;
extra21=zeros(1,n);
extra22=zeros(1,n);
extra23=zeros(1,n);
for j=1:n
    extra21(1,j)=sqrt(x_filter(1,2,j)-data_measurement1(2,1,j))^2+(x_filter(3,2,j)-data_measurement1(2,3,j))^2;
    extra22(1,j)=sqrt((data_measurement(2,1,j)-data_measurement1(2,1,j))^2+(data_measurement(2,2,j)-data_measurement1(2,3,j))^2);
    extra23(1,j)=extra22(1,j)/extra21(1,j);
end
plot(1:n,extra23(1,:),'k:'); 
xlabel('times'),ylabel('RMSE of b');
grid;
end




function [U] = JPDA_func(Q2,m1,ellipse_Volume,y,Z_predic,S)
    global Pd c;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  4.产生互联矩阵A_matrix,其中num表示可行联合事件个数  %%%
    %%%    即对确认矩阵进行拆分，得到互联矩阵，互联矩阵的特性见《雷达数据处理及应用》P130
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_matrix=zeros(m1,3,10000);
        A_matrix(:,1,1:10000)=1;
        if m1~=0      %m1=0表示两个目标都没有观测
           num=1;
           %对所有落入门限内的目标逐一循环判断
           for i=1:m1
                %为1表示该量测点落入目标1门限内
                if Q2(i,2)==1
                    A_matrix(i,2,num)=1;  %属于目标1
                    A_matrix(i,1,num)=0;  %不是虚假目标
                    num=num+1;
                    for j=1:m1
                        if (i~=j)&&(Q2(j,3)==1)
                            %量测点属于目标1
                            A_matrix(i,2,num)=1;
                            A_matrix(i,1,num)=0;

                            %量测点属于目标2
                            A_matrix(j,3,num)=1;
                            A_matrix(j,1,num)=0;
                            num=num+1;
                        end
                    end
                end
           end                                   

            %程序运行到这时，表明目标1存在关联量测的情况已经全部讨论完毕
            for i=1:m1
                if Q2(i,3)==1
                    A_matrix(i,3,num)=1; %属于目标2
                    A_matrix(i,1,num)=0; %不是虚假目标
                    num=num+1;
                end
            end
        else
            flag=1;
        end

        %穷举法拆分的结果存在A_matrix中
        A_matrix=A_matrix(:,:,1:num); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  5.计算后验概率Pr,其中False_num表示假量测,mea_indicator表示观测指示器,target_indicator表示目标指示器
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % m1是确认过与当前目标互联的量测点数量
    % num为可行互联矩阵个数
    Pr=zeros(1,num);
    for i=1:num
        False_num=m1;
        N=1;

        %% 下列步骤计算《雷达数据。。。》P161（8.112）
        for j=1:m1
            %量测互联指示
            mea_indicator=sum(A_matrix(j,2:3,i));  
            if mea_indicator==1  %表明该量测点成功与某一目标互联
                False_num=False_num-1;
                if A_matrix(j,2,i)==1  
                    %如果观测与目标1关联
                    %计算公式见《雷达数据。。。》P125（7.35）
                    b=(y(:,j)-Z_predic(:,1))'*inv(S(:,:,1))*(y(:,j)-Z_predic(:,1));
                    N=N/sqrt(det(2*pi*S(:,:,1)))*exp(-1/2*b);                          %计算正态分布函数                         
                else
                    %如果观测与目标2关联
                    b=(y(:,j)-Z_predic(:,2))'*inv(S(:,:,2))*(y(:,j)-Z_predic(:,2));
                    N=N/sqrt(det(2*pi*S(:,:,2)))*exp(-1/2*b);                          %计算正态分布函数                         
                end                                                                        
            end
        end

        %表示整个空域的体积
        V=ellipse_Volume(1)+ellipse_Volume(2);  
        Pr(i)=N/(V^False_num);

        %% 下列步骤计算《雷达数据。。。》P161（8.114）
        if Pd==1
            a=1;
        else
            a=1;
            for j=1:c
                target_indicator=sum(A_matrix(:,j+1,i));                               %参考文献中式4-49
                %详见《雷达数据。。。》P163(8.116)前半部分
                a=a*Pd^target_indicator*(1-Pd)^(1-target_indicator);                   %计算检测概率
            end
        end  

        a1=1;
        for j=1:False_num
            a1=a1*j;
        end

        %计算每一个量测目标的后验概率

        %详见《雷达数据。。。》P134 公式7.82
        Pr(i)=Pr(i)*a*a1;
    end

    %% 
    %详见《雷达数据。。。》P162(8.109)
    Pr=Pr/sum(Pr);

    %% 6.计算全部确认矩阵和各个目标的关联概率U  
    U=zeros(m1+1,c);
    for i=1:c
        for j=1:m1%量测个数
            for k=1:num %互联矩阵个数（联合事件个数）
                %《雷达数据。。。》P159(8.85)
                U(j,i)=U(j,i)+Pr(k)*A_matrix(j,i+1,k);
            end
        end
    end

    %无量测与目标T互联的关联概率存入U（m1+1,:),归一化
    U(m1+1,:)=1-sum(U(1:m1,1:c)); 
end

  

