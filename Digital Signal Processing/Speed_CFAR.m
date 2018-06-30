function [DotNum,DotLocal] = Speed_CFAR(L_r,AmbData)
DotNum = 0;
N_guard = 2;
N_ref = 6;
L_v = 64;

DotLocal(1) = 0;
DotLocal(2) = 0;
OS_CFAR_Factor = 6.5;
% figure(3)
% clf
for i=1:L_r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OS_CFAR�㷨%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp=AmbData(i,:);
    %��1�㵽��N_guard+1����龯���������ɺ���Nref�������
    for r_n = 1:N_guard+1
        Os_Index = 0;
        for j=1+N_guard+1:N_guard+N_ref
            Os_Index = Os_Index + 1;
            temp_OS_1(Os_Index) = temp(j);
        end

        for kk=1:Os_Index
            [max_os,index_os] = max(temp_OS_1);
            temp_OS_1(index_os) = 0;
            temp_OS_2(Os_Index-kk+1) = max_os;
        end

        %ȥ��OS��������������ֵ
        temp_OS_2(Os_Index) = 0;
        temp_OS_2(Os_Index-1) = 0;
        os_mean = sum(temp_OS_2)/(Os_Index-2);
        OS_CFAR_gate(r_n) = os_mean*OS_CFAR_Factor;   

        if temp(r_n) > OS_CFAR_gate(r_n)
            DotLocal(DotNum*2+1) = i;  %%���뵥Ԫ
            DotLocal(DotNum*2+2) = r_n; %%�ٶȵ�Ԫ
            DotNum = DotNum + 1;
%             plot(r_n,temp(r_n),'ro');
%             hold on ; 
        end
        temp_OS_1 = zeros(1,16);
        temp_OS_2 = zeros(1,16);
    end

    %��Nguard+1�㵽��Nref+Nguard-1����龯����ʱ������ֵ��ǰ���ref���㹲ͬ����
    for r_n = N_guard+2:N_guard+N_ref
        Os_Index = 0;
        for j=1:r_n-N_guard-1
            Os_Index = Os_Index + 1;
            temp_OS_1(Os_Index) = temp(j);         
        end
        for j=r_n+N_guard+1:r_n+N_guard+N_ref
            Os_Index = Os_Index + 1;
            temp_OS_1(Os_Index) = temp(j);
        end

        for kk=1:Os_Index
            [max_os,index_os] = max(temp_OS_1);
            temp_OS_1(index_os) = 0;
            temp_OS_2(Os_Index-kk+1) = max_os;
        end

        temp_OS_2(Os_Index) = 0;
        temp_OS_2(Os_Index-1) = 0;
        os_mean = sum(temp_OS_2)/(Os_Index-2);

        OS_CFAR_gate(r_n) = os_mean*OS_CFAR_Factor; 
        
        if temp(r_n) > OS_CFAR_gate(r_n)
            DotLocal(DotNum*2+1) = i;  %%���뵥Ԫ
            DotLocal(DotNum*2+2) = r_n; %%�ٶȵ�Ԫ
            DotNum = DotNum + 1;
%             plot(r_n,temp(r_n),'ro');
%             hold on ; 
        end
        temp_OS_1 = zeros(1,16);
        temp_OS_2 = zeros(1,16);
    end

    %������Ĵ���
    for r_n=N_ref+N_guard+1:L_v-N_ref-N_guard
        Os_Index = 0;   
        for j=r_n-(N_ref+N_guard):r_n-N_guard-1
            Os_Index = Os_Index + 1;
            temp_OS_1(Os_Index) = temp(j);
        end

        for j=r_n+N_guard+1:r_n+(N_ref+N_guard)
            Os_Index = Os_Index + 1;
            temp_OS_1(Os_Index) = temp(j);
        end

        for kk=1:Os_Index
            [max_os,index_os] = max(temp_OS_1);
            temp_OS_1(index_os) = 0;
            temp_OS_2(Os_Index-kk+1) = max_os;
        end

        temp_OS_2(Os_Index) = 0;
        temp_OS_2(Os_Index-1) = 0;
        os_mean = sum(temp_OS_2)/(Os_Index-2);

        OS_CFAR_gate(r_n) = os_mean*OS_CFAR_Factor; 
        
        if temp(r_n) > OS_CFAR_gate(r_n)
            DotLocal(DotNum*2+1) = i;  %%���뵥Ԫ
            DotLocal(DotNum*2+2) = r_n; %%�ٶȵ�Ԫ
            DotNum = DotNum + 1;
%             plot(r_n,temp(r_n),'ro');
%             hold on ; 
        end
        temp_OS_1 = zeros(1,16);
        temp_OS_2 = zeros(1,16);
    end

    %��Ns-Nguard-Nref�㵽��Ns-1-Nguard-1���������ǰ��Nref���㹲ͬ����
    for r_n = L_v-N_guard-N_ref+1:L_v-N_guard-1
        Os_Index = 0;
        
        for j=r_n-(N_ref+N_guard):r_n-N_guard-1
            Os_Index = Os_Index + 1;
            temp_OS_1(Os_Index) = temp(j);
        end    
        
        for j = r_n+N_guard+1:L_v
            Os_Index = Os_Index + 1;
            temp_OS_1(Os_Index) = temp(j);
        end  
        
        for kk=1:Os_Index
            [max_os,index_os] = max(temp_OS_1);
            temp_OS_1(index_os) = 0;
            temp_OS_2(Os_Index-kk+1) = max_os;
        end
        
        temp_OS_2(Os_Index) = 0;
        temp_OS_2(Os_Index-1) = 0;
        os_mean = sum(temp_OS_2)/(Os_Index-2);
        OS_CFAR_gate(r_n) = os_mean*OS_CFAR_Factor; 
        
        if temp(r_n) > OS_CFAR_gate(r_n)
            DotLocal(DotNum*2+1) = i;  %%���뵥Ԫ
            DotLocal(DotNum*2+2) = r_n; %%�ٶȵ�Ԫ
            DotNum = DotNum + 1;
%             plot(r_n,temp(r_n),'ro');
%             hold on ; 
        end
        temp_OS_1 = zeros(1,16);
        temp_OS_2 = zeros(1,16);
    end

    %��Ns-Nguard-1�㵽��Ns-1��������ֵ����ǰ���Nref�������
    for r_n = L_v-N_guard:L_v
        Os_Index = 0;
        
        for j=r_n-(N_ref*2+N_guard):r_n-N_guard-1
            Os_Index = Os_Index + 1;
            temp_OS_1(Os_Index) = temp(j);
        end

        for kk=1:Os_Index
            [max_os,index_os] = max(temp_OS_1);
            temp_OS_1(index_os) = 0;
            temp_OS_2(Os_Index-kk+1) = max_os;
        end
        
        temp_OS_2(Os_Index) = 0;
        temp_OS_2(Os_Index-1) = 0;
        os_mean = sum(temp_OS_2)/(Os_Index-2);
        OS_CFAR_gate(r_n) = os_mean*OS_CFAR_Factor; 

        if temp(r_n) > OS_CFAR_gate(r_n)
            DotLocal(DotNum*2+1) = i;  %%���뵥Ԫ
            DotLocal(DotNum*2+2) = r_n; %%�ٶȵ�Ԫ
            DotNum = DotNum + 1;
%             plot(r_n,temp(r_n),'ro');
%             hold on ; 
        end
        
        temp_OS_1 = zeros(1,16);
        temp_OS_2 = zeros(1,16);
    end 
    
%     plot(1:64,temp,'-b');
%     hold on;
%     plot(1:64,OS_CFAR_gate,'-.r');
%     hold on 

end
end
            
    