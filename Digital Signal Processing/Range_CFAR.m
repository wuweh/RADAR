function [targetNum,local] = Range_CFAR(dotnum,DotLocal,data)
L_r = size(data,2);
N_guard = 2;
N_ref = 8;
flag = 0;
targetNum = 0;
local(1) = 0;
OS_CFAR_Factor = 4.5;
% figure(2)
% clf
for i=1:dotnum
    rangeIndex = DotLocal((i-1)*2+1);
    speedIndex = DotLocal((i-1)*2+2);
    flag = 0;
    temp = data(speedIndex,:);
    OS_CFAR_gate = zeros(1,64);
    
    %��1�㵽��N_guard+1����龯���������ɺ���Nref�������
    if rangeIndex <= N_guard+1
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
        OS_CFAR_gate(rangeIndex) = os_mean*OS_CFAR_Factor;   

        if temp(rangeIndex) > OS_CFAR_gate(rangeIndex)
            flag = 1;
        end
        temp_OS_1 = zeros(1,16);
        temp_OS_2 = zeros(1,16);
    end
    
    if ( rangeIndex>N_guard+1) && (rangeIndex<=N_ref+N_guard)
        Os_Index = 0;
        for j=1:rangeIndex-N_guard-1
            Os_Index = Os_Index + 1;
            temp_OS_1(Os_Index) = temp(j);
        end
        for j=rangeIndex+N_guard+1:rangeIndex+N_guard+N_ref
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
        OS_CFAR_gate(rangeIndex) = os_mean*OS_CFAR_Factor;   

        if temp(rangeIndex) > OS_CFAR_gate(rangeIndex)
            flag = 1;
        end
        temp_OS_1 = zeros(1,16);
        temp_OS_2 = zeros(1,16);
    end
    
    %������Ĵ���
    if (rangeIndex<=(L_r-N_ref-N_guard)) && (rangeIndex>=(N_ref+N_guard+1))
        Os_Index = 0;
        for j=rangeIndex-(N_ref+N_guard):rangeIndex-N_guard-1
            Os_Index = Os_Index + 1;
            temp_OS_1(Os_Index) = temp(j);
        end

        for j=rangeIndex+N_guard+1:rangeIndex+(N_ref+N_guard)
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
        OS_CFAR_gate(rangeIndex) = os_mean*OS_CFAR_Factor;   

        if temp(rangeIndex) > OS_CFAR_gate(rangeIndex)
            flag = 1;
        end
        temp_OS_1 = zeros(1,16);
        temp_OS_2 = zeros(1,16);
    end
    
    %����ά�����
    if(1==flag)
       local(targetNum*2+1) = speedIndex;
       local(targetNum*2+2) = rangeIndex;
       targetNum = targetNum + 1;     
    end
    
    if targetNum == 0
        local(1) = 0;
        local(2) = 0;
    end
    
%     plot(1:64,temp,'-b');
%     hold on;
%     plot(1:64,OS_CFAR_gate,'-.r');
%     hold on 
%     
%     %���������
%     if temp(rangeIndex)>OS_CFAR_gate(rangeIndex)
%         plot(rangeIndex,temp(rangeIndex),'ro');
%         hold on ; 
%     end  
end
end

    