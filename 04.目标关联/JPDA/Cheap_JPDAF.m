clc;
clear all;
close all;

global Pd Pg g_sigma lambda c gamma target_position_real;
Pd=1;         
Pg=0.97;      
g_sigma=9.21; 
lambda=8;     
c=3;                        %target numbers       
gamma=lambda*10^(-6); 

n=50;   
T=1;    
MC_number=5;     %Monte Carlo
                                    
A = [1 T 0 0;
        0 1 0 0;
        0 0 1 T;
        0 0 0 1];  
    
target_position=[1500 300 800 400; 
                          1000 400 3000 300;
                          9000 100 2000 400];    
target_position_real = zeros(3,4,50);
target_position_real(:,:,1)=target_position;  

%生成航迹
for i=2:n
    for j=1:3
        target_position_real(j,:,i) =  A*target_position_real(j,:,i-1)';
    end
end

% x(1:50)= target_position_real(1,1,1:50);
% y(1:50)= target_position_real(1,3,1:50);
% plot(x,y,'-.');hold on;
% x(1:50)= target_position_real(2,1,1:50);
% y(1:50)= target_position_real(2,3,1:50);
% plot(x,y,'-.');hold on;
% x(1:50)= target_position_real(3,1,1:50);
% y(1:50)= target_position_real(3,3,1:50);
% plot(x,y,'-.');hold on;

JPDAF(target_position,n,T,MC_number,c);     

function  JPDAF(target_position,n,T,MC_number,c)
    global g_sigma gamma target_position_real;                                                                                                             
    target_delta=[100 100 100];                                                                               
    P=zeros(4,4,c);                                                                   
    P1=[target_delta(1)^2 0 0 0;0 0.01 0 0;0 0 target_delta(1)^2 0;0 0 0 0.01];       
    P(:,:,1)=P1;
    P(:,:,2)=P1;
    A = [1 T 0 0;
            0 1 0 0;
            0 0 1 T;
            0 0 0 1];   
    C = [1 0 0 0;
            0 0 1 0];  

    R=[target_delta(1)^2 0;0 target_delta(1)^2];                                       
    Q=[4 0;0 4];                                                                       
    G=[T^2/2 0;T 0;0 T^2/2;0 T];                                                      
    x_filter=zeros(4,c,n);                                                      
    x_filter1=zeros(4,c,n,MC_number);                                                  
    data_measurement=zeros(c,2,n);                                                    

for M=1:MC_number
    Noise=[];
    for i=1:n
        for j=1:c                                                                      
            data_measurement(j,:,i)= C*target_position_real(j,:,i)'+sqrtm(Q)*randn(2,1)*100;
        end
    end
    S=zeros(2,2,c);
    Z_predic=zeros(2,2);                                                               
    x_predic=zeros(4,2);                                                               
    ellipse_Volume=zeros(1,2);

    NOISE_sum_a=[]; 
    NOISE_sum_b=[];
    NOISE_sum_c=[];

    for t=1:n
        y1=[];
        y=[];
        Noise=[];
        NOISE=[];

        for i=1:c      
            if t~=1
                x_predic(:,i) = A*x_filter(:,i,t-1);      
            else
                x_predic(:,i)=target_position(i,:)';                                                   
            end

            P_predic=A*P(:,:,i)*A'+G*Q*G';                                                
            Z_predic(:,i)=C*x_predic(:,i);                                                 
            R=[target_delta(i)^2 0; 0 target_delta(i)^2];
            S(:,:,i)=C*P_predic*C'+R;   
            ellipse_Volume(i)= pi*g_sigma*sqrt(det(S(:,:,i)));         
             number_returns=floor(ellipse_Volume(i)*gamma+1);        
           % number_returns = 2;       
            side=sqrt((ellipse_Volume(i)*gamma+1)/gamma)/2;            
            Noise_x=data_measurement(i,1,t)+side-2*rand(1,number_returns)*side;  
            Noise_y=data_measurement(i,2,t)+side-2*rand(1,number_returns)*side;    
            Noise=[Noise_x;Noise_y];
            NOISE=[NOISE Noise];

            if i==1
                NOISE_sum_a =[NOISE_sum_a Noise]; 
            end

           if i==2
                NOISE_sum_b =[NOISE_sum_b Noise]; 
           end

           if i==3
                NOISE_sum_c =[NOISE_sum_c Noise]; 
            end
        end

        b=zeros(1,2);
        b(1)=data_measurement(1,1,t);
        b(2)=data_measurement(1,2,t);
        y1=[NOISE b'];               
        b(1)=data_measurement(2,1,t);
        b(2)=data_measurement(2,2,t);
        y1=[y1 b'];                  
        b(1)=data_measurement(3,1,t);
        b(2)=data_measurement(3,2,t);
        y1=[y1 b']; 
        [n1,n2]=size(y1);
        
        m1=0;                       
        [n1,n2]=size(y1);
        Q1=zeros(100,3);
        for j=1:n2 
            flag=0;
            for i=1:c
                d=y1(:,j)-Z_predic(:,i);
                D=d'*inv(S(:,:,i))*d; 
                if D<=g_sigma                                                    
                   flag=1;
                   Q1(m1+1,1)=1;
                   Q1(m1+1,i+1)=1;
                end
            end
            if flag==1   
               y=[y y1(:,j)];                                                     
               m1=m1+1                                                   
            end
        end 

     [x_filter(:,:,t), P] = cheap_JPDA_xxx(x_predic, P ,Z_predic ,S ,3, y, m1, A, C, G*Q*G');
     x_filter1(:,:,t,M)  = x_filter(:,:,t);
     
    end
end

    %%下面是作图部分
    x_filter=sum(x_filter1,4)/MC_number;  
%     figure;
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
        if i ==3
           plot(a(:),b(:),'g*')
        end
        hold on;
    end
    plot(NOISE_sum_a(1,:),NOISE_sum_a(2,:),'r.');
    plot(NOISE_sum_b(1,:),NOISE_sum_b(2,:),'b.');
    plot(NOISE_sum_c(1,:),NOISE_sum_c(2,:),'g.');

    hold on;
    a(1:n) = x_filter(1,1,1:n);
    b(1:n) = x_filter(3,1,1:n);
    plot(a,b,'r--.');
    a(1:n) = x_filter(1,2,1:n);
    b(1:n) = x_filter(3,2,1:n);
    plot(a,b,'b--.');
    a(1:n) = x_filter(1,3,1:n);
    b(1:n) = x_filter(3,3,1:n);
    plot(a,b,'g--.');
    legend('Measure A','Measure B','Measure C','Noise A','Noise B','Noise C','Filter A','Filter B','Filter C');grid;

    figure;
    a=0;
    c1=zeros(c,n);
    for j=1:n
        for i=1:MC_number   
            a=(x_filter1(1,1,j,i)-data_measurement(1,1,j))^2+(x_filter1(3,1,j,i)-data_measurement(1,2,j))^2;
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
    [c_extra,pos]=sort(b_extra);                                                       
    a_extra(2,:)=a_extra(2,pos);                                                      
    a_extra(1,:)=c_extra;                                                              
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
        for i=1:MC_number                                                             
            a=(x_filter1(1,2,j,i)-data_measurement(2,1,j))^2+(x_filter1(3,2,j,i)-data_measurement(2,2,j))^2;
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
    [c_extra,pos]=sort(b_extra);                                                    
    a_extra(2,:)=a_extra(2,pos);                                                       
    a_extra(1,:)=c_extra;                                                             
    str1=num2str(a_extra(2,n));
    str2=num2str(a_extra(1,n));
    str=strcat('\itN=',str1,'\itError=',str2,'(m)');
    text(a_extra(2,n),0.8*a_extra(1,n),str);
    hold on;
    plot([a_extra(2,n) a_extra(2,n)],[0 a_extra(1,n)],'b');
    hold on;
    plot(1:n,c1(2,:),'b:');

    a=0;
    for j=1:n
        for i=1:MC_number                                                      
            a=(x_filter1(1,3,j,i)-data_measurement(3,1,j))^2+(x_filter1(3,3,j,i)-data_measurement(3,2,j))^2;
            c1(3,j)=c1(3,j)+a;
        end
            c1(3,j)=sqrt(c1(3,j)/MC_number);
    end
    temp=c1(3,:);
    a_extra=zeros(2,n);
    b_extra=zeros(1,n);
    c_extra=zeros(1,n);
    a_extra(1,:)=temp;
    a_extra(2,:)=1:1:n;
    b_extra=a_extra(1,:);
    [c_extra,pos]=sort(b_extra);                                                    
    a_extra(2,:)=a_extra(2,pos);                                                      
    a_extra(1,:)=c_extra;                                                             
    str1=num2str(a_extra(2,n));
    str2=num2str(a_extra(1,n));
    str=strcat('\itN=',str1,'\itError=',str2,'(m)');
    text(a_extra(2,n),0.8*a_extra(1,n),str);
    hold on;
    plot([a_extra(2,n) a_extra(2,n)],[0 a_extra(1,n)],'g');
    hold on;
    plot(1:n,c1(3,:),'g:');

    xlabel('times'),ylabel('Mean Square Error/m');
    legend('MAX A','MSE A','MAX B','MSE B','MAX C','MSE C');grid;
end


function [U] = cheap_JPDA(m1,y,Z_predic,S)
    global c
    Gjt = zeros(m1,c);
    for i = 1:c
        temp = 0;
        for j = 1:m1
            v1 = y(:,j)-Z_predic(:,i);
            v = (-0.5*v1'*inv(S(:,:,i))*v1);
            Gjt(j,i) = exp(v)/(2*3.14*sqrt(det(S(:,:,i))));
            temp = temp + Gjt(j,i);
        end
        St(i) =  temp;
    end
    
    B = 0;
    for i = 1:c
        ST = St(i);
        for j = 1:m1
            GJT = Gjt(j,i);
            Sj = 0;
            for k =1:c
                Sj = Sj + Gjt(j,k);
            end        
            U(j,i) = GJT/((ST+Sj-GJT+B));
        end
    end
      
    for i = 1:c
        U(j+1,i) = 1 - sum(U(1:j,i));
    end
end


%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%

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
    
%         for i=1:c                                                                  
%         P_predic = A*P(:,:,i)*A'+G*Q*G';
%         K(:,:,i)= P_predic*C'*inv(S(:,:,i));
%         P(:,:,i)= P_predic-(1-U(m1+1,i))*K(:,:,i)*S(:,:,i)*K(:,:,i)';
%     end
% 
%     for i=1:c
%         a=0;         
%         b=0;
% 
%         x_filter2=0;
%         for j=1:m1 
%             x_filter2=x_filter2+U(j,i)*(x_predic(:,i)+ K (:,:,i)*(y(:,j)- Z_predic(:,i)));
%         end
%         x_filter2=U(j+1,i)*x_predic(:,i)+x_filter2;
%         x_filter(:,i,t)=x_filter2;
% 
%         for j=1:m1+1
%             if j==m1+1
%                 a=x_predic(:,i);
%             else
%                a=x_predic(:,i)+ K (:,:,i)*(y(:,j)- Z_predic(:,i));
%             end
%             b=b+U(j,i)*(a*a'-x_filter2*x_filter2');
%         end
%         P(:,:,i)=P(:,:,i)+b; 
%         x_filter1(:,i,t,M)=x_filter(:,i,t);
%     end
    
end

