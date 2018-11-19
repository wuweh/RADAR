function [x_putput, P] = PDA_Function(x_predic, P_predic, S, Z_predic,K, y, m)
    global Pd Pg gamma hfun;
    % 计算各个杂波的权重
    % gamma为空间杂波密度，在实际工程中
    Bk=gamma*sqrt(det(S)*2*3.14)*(1-Pd*Pg)/Pd; %算b     
    %m=0表示无有效回波 
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
        belta0=Bk/(Bk+sum(E));    %计算b0  %无回波时的关联概率 
        v=zeros(4,1); 
        v1=zeros(4,1); 
        for i=1:m 
            belta(i)=E(i)/(Bk+sum(E));     %算关联概率 计算bi
            v=v+belta(i)*(y(:,i)-Z_predic); 
            v1=v1+belta(i)*(y(:,i)-Z_predic)*(y(:,i)-Z_predic)'; 
        end 
        %计算预测状态值
        x_filter(:,1)= x_predic + K*v; 
        %计算预测状态协方差 
        Pc=(eye(4)-K*hfun)*P_predic; 
        PP=K*(v1-v*v')*K'; 
        P=belta0*P_predic+(1-belta0)*Pc+PP;
    end
    x_putput = x_filter(:,1);
end
