clc;
clear all;
close all;

global Pd Pg g_sigma lambda c gamma;
Pd=1;         %������
Pg=0.99;      %��ȷ��������������ڵø���
g_sigma=9.21; %���ޣ�����ֵ��
lambda=8;     %�Ӳ��ܶ�
c=2;          %Ŀ�����
gamma=lambda*10^(-6); 

n=50;   %��������
T=1;    %TΪ�������
MC_number=5;     %Monte Carlo�������
                                                                             
target_position=[1500 300 800 400; 1000 400 3000 300];        %Ŀ�����ʼλ�ú��ٶ�(m,m/s)                   
JPDAF(target_position,n,T,MC_number,c);     


function  JPDAF(target_position,n,T,MC_number,c)
global g_sigma gamma;                                                             %ÿһ����λ���(km^2)�ڲ���lambda���Ӳ�
% Target_measurement=zeros(c,2,n);                                                   %Ŀ��۲⻥���洢����
target_delta=[100 100 100];                                                            %Ŀ���Ӧ�Ĺ۲��׼��                    
P=zeros(4,4,c);                                                                    %Э�������
P1=[target_delta(1)^2 0 0 0;0 0.01 0 0;0 0 target_delta(1)^2 0;0 0 0 0.01];        %��ʼЭ������� 
P(:,:,1)=P1;
P(:,:,2)=P1;
A = [1 T 0 0;
    0 1 0 0;
    0 0 1 T;
    0 0 0 1];   %״̬ת�ƾ���
C = [1 0 0 0;
    0 0 1 0];  %�۲����

R=[target_delta(1)^2 0;0 target_delta(1)^2];                                       %�۲�Э�������
Q=[4 0;0 4];                                                                       %ϵͳ��������Э����
G=[T^2/2 0;T 0;0 T^2/2;0 T];                                                       %������������
x_filter=zeros(4,c,n);                                                             %�洢Ŀ��ĸ�ʱ�̵��˲�ֵ
x_filter1=zeros(4,c,n,MC_number);                                                  %MC_number��Montle Carlo��������ȫ������洢����
data_measurement=zeros(c,2,n);                                                     %�۲�洢����
data_measurement1=zeros(c,4,n);                                                    %ʵ��λ������x,y����   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  ����Ŀ���ʵ��λ��  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_measurement1(:,:,1)=target_position;                                          %ʵ��λ�þ����ʼ�� 
for i=1:c
    for ii=2:n    
        %��������Ŀ���ʵ��λ�� 
        data_measurement1(i,:,ii)=(A*data_measurement1(i,:,ii-1)')'+(G*sqrt(Q)*(randn(2,1)))';        
    end
end

for M=1:MC_number
%%%%%%%%%%%%%%%%%%%%    
%%%  1.����·��  %%%
%%%%%%%%%%%%%%%%%%%%
Noise=[];
for i=1:n
    for j=1:c                                                                      
        %�õ�����Ŀ��Ĺ۲�λ�ã���������
        data_measurement(j,1,i)=data_measurement1(j,1,i)+rand(1)*target_delta(j);
        data_measurement(j,2,i)=data_measurement1(j,3,i)+rand(1)*target_delta(j); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  2.�����Ӳ�,��ȷ����Ч�۲�  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=zeros(2,2,c);

%�洢����Ŀ��Ĺ۲�Ԥ��ֵ,��ֻ����x,y����
Z_predic=zeros(2,2);                                                               

%�洢����Ŀ���״̬Ԥ��ֵ,������x,y�����x,y�����ٶ�
x_predic=zeros(4,2);                                                               

ellipse_Volume=zeros(1,2);

 %�洢Ŀ��1���Ӳ�
NOISE_sum_a=[]; 
 %�洢Ŀ��2���Ӳ�
NOISE_sum_b=[];
 %�洢Ŀ��3���Ӳ�
NOISE_sum_c=[];

for t=1:n
    y1=[];
    y=[];
    Noise=[];
    NOISE=[];
    
    for i=1:c      
        if t~=1
            %��ǰһʱ�̵��˲�ֵ��Ԥ�⵱ǰ��ֵ(kalman�˲��ĵ�һ�����ʽ)
            x_predic(:,i) = A*x_filter(:,i,t-1);      
%             x_predic(:,i)=data_measurement1(i,:,t)'; 
        else
            %��һ�β�����������ʵλ�õ�Ԥ��ֵ
            x_predic(:,i)=target_position(i,:)';                                                   
        end
        
        %����x_predicЭ�������(kalman�˲��ĵڶ������ʽ) 
        P_predic=A*P(:,:,i)*A'+G*Q*G';                                                
        
        %��ȡԤ��ֵ��x,y���꣬����x,y�ٶ�                                  
        Z_predic(:,i)=C*x_predic(:,i);                                                 
        R=[target_delta(i)^2 0; 0 target_delta(i)^2];
    
        %SΪ��ϢЭ����
        S(:,:,i)=C*P_predic*C'+R;   
        
        %���в�����Ϊ�����ɶ���ز����������Ϊ���Ŀ��
        %����ȷ���������
        ellipse_Volume(i)= pi*g_sigma*sqrt(det(S(:,:,i)));          %������Բ�����ŵ����   
%          number_returns=floor(ellipse_Volume(i)*gamma+1);           %��Բ�������ڵĴ���ز���
        number_returns = 2;       
        side=sqrt((ellipse_Volume(i)*gamma+1)/gamma)/2;            %����Բ�����ŵ�ЧΪ�����Σ�����������α߳��Ķ���֮һ
        Noise_x=x_predic(1,i)+side-2*rand(1,number_returns)*side;  %��Ԥ��ֵ��Χ��������ز���ע�⣺��ĳһ��number_returnsС�ڵ���0ʱ�����������һ�μ��ɡ�
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
    y1=[NOISE b'];                %�����յ������еĻز�����y1��,�����Ӳ��͹۲�
    b(1)=data_measurement(2,1,t);
    b(2)=data_measurement(2,2,t);
    y1=[y1 b'];                  
    [n1,n2]=size(y1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  3.�����۲�ȷ�Ͼ���Q2  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m1=0;                         %��¼��Ч�۲����
    [n1,n2]=size(y1);
    Q1=zeros(100,3);
    for j=1:n2 
        flag=0;
        for i=1:c
            %��������ֵ��Ԥ��ֵ�Ĳв�
            d=y1(:,j)-Z_predic(:,i);
            D=d'*inv(S(:,:,i))*d; 
            %�ж��Ƿ���������Ҫ��������ʵ�ǰ��տ����ֲ�������ʣ�
            if D<=g_sigma                                                    
               flag=1;
               Q1(m1+1,1)=1;
               Q1(m1+1,i+1)=1;
            end
        end
        if flag==1   
           y=[y y1(:,j)];                                                      %������������е����лز�����y��
           m1=m1+1;                                                            %��¼��Ч�۲����
        end
    end
    Q2=Q1(1:m1,1:3);
    [U] = JPDA_func(Q2,m1,ellipse_Volume,y,Z_predic,S);

%%  7.Kalman�˲���ʼ  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:c                                                                          %����Э�������
    P_predic = A*P(:,:,i)*A'+G*Q*G';
    %����������
    K(:,:,i)= P_predic*C'*inv(S(:,:,i));
    %������״����ݡ�������P165(8.130����1)
    P(:,:,i)= P_predic-(1-U(m1+1,i))*K(:,:,i)*S(:,:,i)*K(:,:,i)';
end

%% 
for i=1:c
    a=0;         
    b=0;
    
    %% ����״̬Ԥ��ֵ
    x_filter2=0;
    for j=1:m1 
        x_filter2=x_filter2+U(j,i)*(x_predic(:,i)+ K (:,:,i)*(y(:,j)- Z_predic(:,i)));
    end
    x_filter2=U(j+1,i)*x_predic(:,i)+x_filter2;
    x_filter(:,i,t)=x_filter2;
    
    %% ����״̬Э����
    for j=1:m1+1
        if j==m1+1
            a=x_predic(:,i);
        else
            %����Ԥ��ֵ
           a=x_predic(:,i)+ K (:,:,i)*(y(:,j)- Z_predic(:,i));
        end
        %������״����ݡ�������P135 ��ʽ7.102�ڶ�����
        b=b+U(j,i)*(a*a'-x_filter2*x_filter2');
    end
    %Э����
    %������״����ݡ�������P136 ��ʽ7.102
    P(:,:,i)=P(:,:,i)+b; 
    x_filter1(:,i,t,M)=x_filter(:,i,t);
end
end
end
 

%% ��ͼ����
% x_filter=sum(x_filter1,4)/MC_number;   %�˲�ֵ��ƽ��
%%%%%%%%%%%%%%%%%%%%
%%%  1.�˲����  %%%
%%%%%%%%%%%%%%%%%%%%
figure;
%Ŀ��a,b�Ĺ۲�λ��
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
%Ŀ��a,b���Ӳ�λ��
plot(NOISE_sum_a(1,:),NOISE_sum_a(2,:),'r.');
plot(NOISE_sum_b(1,:),NOISE_sum_b(2,:),'b.');

hold on;
%Ŀ��a,b�Ĺ���λ��
a(1:n) = x_filter(1,1,1:n);
b(1:n) = x_filter(3,1,1:n);
plot(a,b,'r-');
a(1:n) = x_filter(1,2,1:n);
b(1:n) = x_filter(3,2,1:n);
plot(a,b,'b-');
legend('Ŀ��a�Ĺ۲�λ��','Ŀ��b�Ĺ۲�λ��','Ŀ��a���Ӳ�','Ŀ��b���Ӳ�','Ŀ��a�Ĺ���λ��','Ŀ��b�Ĺ���λ��');grid;


%%%%%%%%%%%%%%%%%%%%
%%%  2.�ٶ����  %%%
%%%%%%%%%%%%%%%%%%%%
figure;
a=0;
c1=zeros(c,n);
for j=1:n
    for i=1:MC_number    %��С�������
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
[c_extra,pos]=sort(b_extra);                                                       %posΪ�������±�,cΪ��һ�е�������;
a_extra(2,:)=a_extra(2,pos);                                                       %�ڶ��а��յ�һ��������±��Ӧ
a_extra(1,:)=c_extra;                                                              %��һ�н�����¸���a �ĵ�һ��;
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
    for i=1:MC_number                                                              %��С�������
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
[c_extra,pos]=sort(b_extra);                                                       %posΪ�������±�,cΪ��һ�е�������;
a_extra(2,:)=a_extra(2,pos);                                                       %�ڶ��а��յ�һ��������±��Ӧ
a_extra(1,:)=c_extra;                                                              %��һ�н�����¸���a �ĵ�һ��;
str1=num2str(a_extra(2,n));
str2=num2str(a_extra(1,n));
str=strcat('\itN=',str1,'\itError=',str2,'(m)');
text(a_extra(2,n),0.8*a_extra(1,n),str);
hold on;
plot([a_extra(2,n) a_extra(2,n)],[0 a_extra(1,n)],'b');
hold on;
plot(1:n,c1(2,:),'b:');

xlabel('times'),ylabel('����ֵ�����ֵ������/m');
legend('Ŀ��a��������ֵ','Ŀ��a�����','Ŀ��b��������ֵ','Ŀ��b�����');grid;

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
    %%%  4.������������A_matrix,����num��ʾ���������¼�����  %%%
    %%%    ����ȷ�Ͼ�����в�֣��õ��������󣬻�����������Լ����״����ݴ���Ӧ�á�P130
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A_matrix=zeros(m1,3,10000);
        A_matrix(:,1,1:10000)=1;
        if m1~=0      %m1=0��ʾ����Ŀ�궼û�й۲�
           num=1;
           %���������������ڵ�Ŀ����һѭ���ж�
           for i=1:m1
                %Ϊ1��ʾ�����������Ŀ��1������
                if Q2(i,2)==1
                    A_matrix(i,2,num)=1;  %����Ŀ��1
                    A_matrix(i,1,num)=0;  %�������Ŀ��
                    num=num+1;
                    for j=1:m1
                        if (i~=j)&&(Q2(j,3)==1)
                            %���������Ŀ��1
                            A_matrix(i,2,num)=1;
                            A_matrix(i,1,num)=0;

                            %���������Ŀ��2
                            A_matrix(j,3,num)=1;
                            A_matrix(j,1,num)=0;
                            num=num+1;
                        end
                    end
                end
           end                                   

            %�������е���ʱ������Ŀ��1���ڹ������������Ѿ�ȫ���������
            for i=1:m1
                if Q2(i,3)==1
                    A_matrix(i,3,num)=1; %����Ŀ��2
                    A_matrix(i,1,num)=0; %�������Ŀ��
                    num=num+1;
                end
            end
        else
            flag=1;
        end

        %��ٷ���ֵĽ������A_matrix��
        A_matrix=A_matrix(:,:,1:num); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  5.����������Pr,����False_num��ʾ������,mea_indicator��ʾ�۲�ָʾ��,target_indicator��ʾĿ��ָʾ��
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % m1��ȷ�Ϲ��뵱ǰĿ�껥�������������
    % numΪ���л����������
    Pr=zeros(1,num);
    for i=1:num
        False_num=m1;
        N=1;

        %% ���в�����㡶�״����ݡ�������P161��8.112��
        for j=1:m1
            %���⻥��ָʾ
            mea_indicator=sum(A_matrix(j,2:3,i));  
            if mea_indicator==1  %�����������ɹ���ĳһĿ�껥��
                False_num=False_num-1;
                if A_matrix(j,2,i)==1  
                    %����۲���Ŀ��1����
                    %���㹫ʽ�����״����ݡ�������P125��7.35��
                    b=(y(:,j)-Z_predic(:,1))'*inv(S(:,:,1))*(y(:,j)-Z_predic(:,1));
                    N=N/sqrt(det(2*pi*S(:,:,1)))*exp(-1/2*b);                          %������̬�ֲ�����                         
                else
                    %����۲���Ŀ��2����
                    b=(y(:,j)-Z_predic(:,2))'*inv(S(:,:,2))*(y(:,j)-Z_predic(:,2));
                    N=N/sqrt(det(2*pi*S(:,:,2)))*exp(-1/2*b);                          %������̬�ֲ�����                         
                end                                                                        
            end
        end

        %��ʾ������������
        V=ellipse_Volume(1)+ellipse_Volume(2);  
        Pr(i)=N/(V^False_num);

        %% ���в�����㡶�״����ݡ�������P161��8.114��
        if Pd==1
            a=1;
        else
            a=1;
            for j=1:c
                target_indicator=sum(A_matrix(:,j+1,i));                               %�ο�������ʽ4-49
                %������״����ݡ�������P163(8.116)ǰ�벿��
                a=a*Pd^target_indicator*(1-Pd)^(1-target_indicator);                   %���������
            end
        end  

        a1=1;
        for j=1:False_num
            a1=a1*j;
        end

        %����ÿһ������Ŀ��ĺ������

        %������״����ݡ�������P134 ��ʽ7.82
        Pr(i)=Pr(i)*a*a1;
    end

    %% 
    %������״����ݡ�������P162(8.109)
    Pr=Pr/sum(Pr);

    %% 6.����ȫ��ȷ�Ͼ���͸���Ŀ��Ĺ�������U  
    U=zeros(m1+1,c);
    for i=1:c
        for j=1:m1%�������
            for k=1:num %������������������¼�������
                %���״����ݡ�������P159(8.85)
                U(j,i)=U(j,i)+Pr(k)*A_matrix(j,i+1,k);
            end
        end
    end

    %��������Ŀ��T�����Ĺ������ʴ���U��m1+1,:),��һ��
    U(m1+1,:)=1-sum(U(1:m1,1:c)); 
end

  

