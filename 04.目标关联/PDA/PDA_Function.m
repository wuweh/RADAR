function [x_putput, P] = PDA_Function(x_predic, P_predic, S, Z_predic,K, y, m)
    global Pd Pg lambda C;
    Bk=lambda*sqrt(det(S)*2*3.14)*(1-Pd*Pg)/Pd;    
    if m==0 
       x_filter(:,1)= x_predic; 
       P=P_predic;    
    else         
        E=zeros(1,m); 
        belta=zeros(1,m); 
        for i=1:m 
            a=(y(:,i)-Z_predic)'*inv(S)*(y(:,i)-Z_predic); 
            E(i)=E(i)+exp(-a/2); 
        end 

        belta0=Bk/(Bk+sum(E));    
        v=zeros(2,1); 
        v1=zeros(2,2); 
        for i=1:m 
            belta(i)=E(i)/(Bk+sum(E));     
            v=v+belta(i)*(y(:,i)-Z_predic); 
            v1=v1+belta(i)*(y(:,i)-Z_predic)*(y(:,i)-Z_predic)'; 
        end 
        x_filter(:,1)= x_predic + K*v; 
        Pc=(eye(4)-K*C)*P_predic; 
        PP=K*(v1-v*v')*K'; 
        P=belta0*P_predic+(1-belta0)*Pc+PP;
    end
    
    x_putput = x_filter(:,1);
end
