function [x_putput, P] = PDA_Function(x_predic, P_predic, S, Z_predic,K, y, m)
    global Pd Pg gamma C;
    % ��������Ӳ���Ȩ��
    %�״����ݴ���Ӧ�� P151(8.42)
    % gammaΪ�ռ��Ӳ��ܶȣ���ʵ�ʹ�����
    Bk=gamma*sqrt(det(S)*2*3.14)*(1-Pd*Pg)/Pd; %��b     
    %m=0��ʾ����Ч�ز� 
    if m==0 
       x_filter(:,1)= x_predic; 
       P=P_predic;    %�޻ز������ 
    else         
        E=zeros(1,m); 
        belta=zeros(1,m); 
        for i=1:m 
            %�״����ݴ���Ӧ�� P151(8.41)
            a=(y(:,i)-Z_predic)'*inv(S)*(y(:,i)-Z_predic); 
            E(i)=E(i)+exp(-a/2); 
        end 
        %�״����ݴ���Ӧ�� P151(8.42)
        belta0=Bk/(Bk+sum(E));    %����b0  %�޻ز�ʱ�Ĺ������� 
        v=zeros(4,1); 
        v1=zeros(4,1); 
        for i=1:m 
            %�״����ݴ���Ӧ�� P151(8.44)
            belta(i)=E(i)/(Bk+sum(E));     %��������� ����bi
            %�״����ݴ���Ӧ�� P148 (8.16)
            v=v+belta(i)*(y(:,i)-Z_predic); 
            %�״����ݴ���Ӧ�� P148��ע1
            v1=v1+belta(i)*(y(:,i)-Z_predic)*(y(:,i)-Z_predic)'; 
        end 
        %�״����ݴ���Ӧ�� P148(8.15)
        %����Ԥ��״ֵ̬
        x_filter(:,1)= x_predic + K*v; 
        %����Ԥ��״̬Э���� 
        %�״����ݴ���Ӧ�� P148(8.18)
        Pc=(eye(4)-K*C)*P_predic; 
        %�״����ݴ���Ӧ�� P148(8.19)
        PP=K*(v1-v*v')*K'; 
        %�״����ݴ���Ӧ�� P148(8.17)
        P=belta0*P_predic+(1-belta0)*Pc+PP;
    end
    
    x_putput = x_filter(:,1);
end
