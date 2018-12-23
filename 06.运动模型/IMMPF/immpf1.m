%状态间交互输入，状态交互输出
function [xe] = immpf1(x,y,R,x0)
T = 1;
N = 100;               %粒子数目
NT = 100;
pai=[0.90,0.10;0.10,0.90]; %定义一步转移概率矩阵
u_cv=0.5;              %匀速运动模型在初始时刻正确的概率
u_ca=0.5;              %匀加速运动模型在初始时刻正确概率
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

for t = 1:NT-1
    
    %1 输入交互
    c_1 = pai(1,1)*u_cv+pai(2,1)*u_ca;
    c_2 = pai(1,2)*u_cv+pai(2,2)*u_ca;
    u11 = pai(1,1)*u_cv/c_1;
    u12 = pai(1,2)*u_cv/c_2;
    u21 = pai(2,1)*u_ca/c_1;
    u22 = pai(2,2)*u_ca/c_2;
    
    x11(:,t) = x1(:,t)*u11+x2(:,t)*u21;
    P11 = u11*(P1+(x1(:,t)-x11(:,t))*(x1(:,t)-x11(:,t))')+u21*(P2+(x2(:,t)-x11(:,t))*(x2(:,t)-x11(:,t))');
    
    x22(:,t) = x1(:,t)*u12+x2(:,t)*u22;
    P22 = u12*(P1+(x1(:,t)-x22(:,t))*(x1(:,t)-x22(:,t))')+u22*(P2+(x2(:,t)-x22(:,t))*(x2(:,t)-x22(:,t))');
    
    %2 产生粒子
    for i = 1:N
        x1part(:,i) = x11(:,t)+sqrtm(P1)*randn(4,1);
        x2part(:,i) = x22(:,t)+sqrtm(P2)*randn(4,1);
    end
    
    %3 PF滤波
    
    %模型1的滤波
    for i = 1:N
        x1partemp(:,i) = A1*x1part(:,i)+15*sqrtm(Q1)*randn(4,1);
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
    
   x1(:,t+1) = [mean(x1part(1,:));mean(x1part(2,:));mean(x1part(3,:));mean(x1part(4,:))];
   y1(:,t+1) = [mean(y1part(1,:));mean(y1part(2,:))];
   P1 = Q1;
   for i = 1:N
       e1(:,i) = x1part(:,i) - x1(:,t+1); 
       P1 = P1 + e1(:,i)*e1(:,i)'*w1(i);
   end
   
    S1 =R;
    for i = 1:N
       S1 = R+w1(i)*er1(:,i)*er1(:,i)';
    end
   
    
    %模型2的滤波
    for i = 1:N
        x2partemp(:,i) = A2*x2part(:,i)+15*sqrtm(Q2)*randn(4,1);
        r2(i) = sqrt(x2partemp(1,i)^2+x2partemp(3,i)^2);
        sita2(i) = atan(x2partemp(3,i)./x2partemp(1,i)); 
        y2part(:,i) = [r2(i);sita2(i)];
        error2(:,i) = y2part(:,i) - y(:,t);
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
    
   x2(:,t+1) = [mean(x2part(1,:));mean(x2part(2,:));mean(x2part(3,:));mean(x2part(4,:))];
   y2(:,t+1) = [mean(y2part(1,:));mean(y2part(2,:))];
    
   P2 = Q2;
   for i = 1:N
       e2(:,i) = x2part(:,i) - x2(:,t+1);
       P2 = P2 + e2(:,i)*e2(:,i)'*w2(i);
   end
   
   S2 =R;
    for i = 1:N
       S2 = R+w2(i)*er2(:,i)*er2(:,i)';
    end                                  
      
   
    %4概率模型更新
        
    Er1 = y(:,t+1) -y1(:,t+1);
    Er2 = y(:,t+1) -y2(:,t+1);
    
    like1=exp(-Er1'*inv(S1)*Er1/2)/(sqrt(2*pi*det(S1)));
    like2=exp(-Er2'*inv(S2)*Er2/2)/(sqrt(2*pi*det(S2)));
    c_1=pai(1,1)*u_cv+pai(2,1)*u_ca;
    c_2=pai(1,2)*u_cv+pai(2,2)*u_ca;
    c=like1*c_1+like2*c_2;
    u_cv=like1*c_1/c; 
    u_ca=like2*c_2/c; 
    
    %5 交互输出
     xe(:,t+1)=x1(:,t+1)*u_cv+x2(:,t+1)*u_ca;
     
end


 