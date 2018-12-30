function [x_putput, P, likehood] = PDA_Function(x_pre, P_pre, fhun, hfun, z_noise)
    global Pd Pg gamma G Q R noise_total;

    x_predic = fhun*x_pre;
    P_predic = fhun*P_pre*fhun'+G*Q*G'; 
    Z_predic = hfun*x_predic;
    S = hfun*P_predic*hfun'+ R; 
    K = P_predic*hfun'*inv(S);
    y=[];
    m=0;
    for j=1:noise_total
        d=z_noise(j,:)'-Z_predic;  
        D(j)=d'*inv(S)*d;  
        if D(j)<=4
            y=[y z_noise(j,:)'];  m=m+1;         
        end 
    end
    z_dem = size(Z_predic,1);
    x_dem = size(x_predic,1);
    v=zeros(z_dem,1); 
    v1=zeros(z_dem,1); 
    Bk=gamma*sqrt(det(S)*2*3.14)*(1-Pd*Pg)/Pd; 
    if m==0 
       x_filter(:,1)= x_predic; 
       P=P_predic;    
    else         
        E=zeros(1,m); 
        belta=zeros(1,m); 
        for i=1:m 
            a=(y(:,i)-Z_predic)'*inv(S)*(y(:,i)-Z_predic); 
            E(i)=exp(-a/2); 
        end 

        belta0=Bk/(Bk+sum(E));    
        for i=1:m 
            belta(i)=E(i)/(Bk+sum(E));     
            v=v+belta(i)*(y(:,i)-Z_predic); 
            v1=v1+belta(i)*(y(:,i)-Z_predic)*(y(:,i)-Z_predic)'; 
        end 
        x_filter(:,1)= x_predic + K*v; 
        Pc=(eye(x_dem)-K*hfun)*P_predic; 
        PP=K*(v1-v*v')*K'; 
        P=belta0*P_predic+(1-belta0)*Pc+PP;
    end
    
    x_putput = x_filter(:,1);
    
    %计算似然函数
    %IMM算法独有
    likehoodmea=[];
    for i=1:m
        likehoodmea(i)=E(i)/(sqrt(2*pi*det(S)));
    end
    temp=[];
    for i=1:m
        LL = likehoodmea(i)*Pg^(-1);
        temp=[temp LL];
    end
    likehood=sum(temp);
end
