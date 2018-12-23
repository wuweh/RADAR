%粒子之间交互且粒子交互输出
function [xe] = immpf2(x,y,R,x0)
T = 1;
N = 100;               %粒子数目
NT = 100;
pai=[0.90,0.10;0.10,0.90]; %定义一步转移概率矩阵
v = 0.21;
P1 =2*diag([1,0.01,1,0.01]);
Q1 =diag([1,0.01,1,0.01]); 
A1 = [1,T,0,0;
      0,1,0,0;
      0,0,1,T;
      0,0,0,1];
  
P2 = 4*diag([2,0.01,2,0.01]);  
Q2 = diag([2,0.01,1,0.01]);
A2 = [1,sin(v*T)/v,0,-(1-cos(v*T))/v;
     0,cos(v*T),0 -sin(v*T);
     0,(1-cos(v*T))/v,1,sin(v*T)/v;
     0,sin(v*T),0,cos(v*T)];
 

x1(:,1) = x0;
x2(:,1) = x0;
xe(:,1) = x0;

%%PARTICLE FILTER 
%%%-----------------------------

for i = 1:N
    x1part(:,i) = x1(:,1)+sqrtm(P1)*randn(4,1);
    x2part(:,i) = x2(:,1)+sqrtm(P2)*randn(4,1);
end
u_cv = ones(1,N)*0.5;
u_ca = ones(1,N)*0.5;
for t = 1:NT-1
    
    %1 输入交互
    for i = 1:N
        c_1(i) = pai(1,1)*u_cv(i)+pai(2,1)*u_ca(i);
        c_2(i) = pai(1,2)*u_cv(i)+pai(2,2)*u_ca(i);
        u11(i) = pai(1,1)*u_cv(i)/c_1(i);
        u12(i) = pai(1,2)*u_cv(i)/c_2(i);
        u21(i) = pai(2,1)*u_ca(i)/c_1(i);
        u22(i) = pai(2,2)*u_ca(i)/c_2(i);
   
        x11part(:,i) = u11(i)*x1part(:,i)+u21(i)*x2part(:,i);
        Pv(:,:,i) = u11(i)*(P1+(x1part(:,i)-x11part(:,i))*(x1part(:,i)-x11part(:,i))')+u21(i)*(P2+(x2part(:,i)-x11part(:,i))*(x2part(:,i)-x11part(:,i))'); 
        
        x22part(:,i) = u12(i)*x1part(:,i)+u22(i)*x2part(:,i);
        Pa(:,:,i) = u12(i)*(P1+(x1part(:,i)-x22part(:,i))*(x1part(:,i)-x22part(:,i))')+u22(i)*(P2+(x2part(:,i)-x22part(:,i))*(x2part(:,i)-x22part(:,i))'); 
    end
    
    %2 滤波
    %模型1的滤波
    for i = 1:N
        x1partemp(:,i) = A1*x11part(:,i)+15*sqrtm(Q1)*randn(4,1);
        r1(i) = sqrt(x1partemp(1,i)^2+x1partemp(3,i)^2);
        sita1(i) = atan(x1partemp(3,i)./x1partemp(1,i)); 
        y1part(:,i) = [r1(i);sita1(i)];
        error1(:,i) = y1part(:,i) - y(:,t+1);
    end
    
    scale11 = max(abs(error1(1,:)))./8;
    scale12 = max(abs(error1(2,:)))./8;
    for i = 1:N
        er1(:,i) = [error1(1,i)/scale11;error1(2,i)/scale12];
        w1(i) = exp(-er1(:,i)'*er1(:,i)/2);
    end
    
    % 权值归一化
    w1sum = sum(w1);
    for i = 1 : N
        w1(i) = w1(i) / w1sum;
    end
               
    c(1) = w1(1);
    for i = 2 : N
        c(i) = c(i-1)+w1(i);
    end
    i = 1;
    u(1) = rand/N;
    for j = 1:N
        u(j) = u(1)+(j-1)/N;
        while u(j)>c(i)
            i = i+1;
        end
        x1part(:,j) = x1partemp(:,i);
        w1(j) = 1/N;
    end
   S2 =zeros(2,2,N);
    for i = 1:N
       S1(:,:,i) = R+er1(:,i)*er1(:,i)';
    end
    
    %模型2的滤波
    for i = 1:N
        x2partemp(:,i) = A2*x2part(:,i)+15*sqrtm(Q2)*randn(4,1);
        r2(i) = sqrt(x2partemp(1,i)^2+x2partemp(3,i)^2);
        sita2(i) = atan(x2partemp(3,i)./x2partemp(1,i)); 
        y2part(:,i) = [r2(i);sita2(i)];
        error2(:,i) = y2part(:,i) - y(:,t+1);
    end
    
    scale21 = max(abs(error2(1,:)))./8;
    scale22 = max(abs(error2(2,:)))./8;
    for i = 1:N
        er2(:,i) = [error2(1,i)/scale21;error2(2,i)/scale22];
        w2(i) = exp(-er2(:,i)'*er2(:,i)/2);
    end
    
    % 权值归一化
    w2sum = sum(w2);
    for i = 1 : N
        w2(i) = w2(i) / w2sum;
    end
     
    c(1) = w2(1);
    for i = 2 : N
        c(i) = c(i-1)+w2(i);
    end
    i = 1;
    u(1) = rand/N;
    for j = 1:N
        u(j) = u(1)+(j-1)/N;
        while u(j)>c(i)
            i = i+1;
        end
        x2part(:,j) = x2partemp(:,i);
        w2(j) = 1/N;
    end
    S2 =zeros(2,2,N);
    for i = 1:N
        S2(:,:,i) =R+er2(:,i)*er2(:,i)';
    end                                                                                                                                                                                                                                                                              
    
    %3 概率模型更新
    for i = 1:N
        
       Er1 = er1(:,i);
       Er2 = er2(:,i);
       
       like1 = exp(-Er1'*inv(S1(:,:,i))*Er1/2)/(sqrt(2*pi));
       like2 = exp(-Er2'*inv(S1(:,:,i))*Er2/2)/(sqrt(2*pi));
       c_1(i)   = pai(1,1)*u_cv(i)+pai(2,1)*u_ca(i);
       c_2(i)   = pai(1,2)*u_cv(i)+pai(2,2)*u_ca(i);
       c(i)     = like1*c_1(i)+like2*c_2(i);
       u_cv(i)  = like1*c_1(i)/c(i); 
       u_ca(i)  = like2*c_2(i)/c(i); 
    end
    
    %5 交互输出
     xe(:,t+1) = zeros(4,1);
     for i = 1:N
         xx(:,i) = w1(i)*x1part(:,i)*u_cv(i)+w2(i)*x2part(:,i)*u_ca(i);
         xe(:,t+1) = xe(:,t+1)+xx(:,i);
     end     
end
