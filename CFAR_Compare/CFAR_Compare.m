


% CFAR算法比较
%**************************************修改记录************************************
%20180624:   1）基于DYC_CFAR算法
%20180624:   1）加入改进后CA_CFAR算法
%20180626:   1）加入OS_CFAR算法
%***************************************END***************************************


clc;
close all;
clear all;

load('ADC_DATA.mat')

DotNum = 0;
N_guard = 2;
N_ref = 6;
L_v = 64;
sum1 = 0;
sum2 = 0;
DotLocal(1) = 0;
DotLocal(2) = 0;

pfa=1e-6;
T=pfa^(-1/N_ref)-1;%均匀杂波背景下标称化因子T
figure;
for i=1:64
    temp=AmbData(i,:);
    avg_whole = mean(temp);
    
    %第1点到第N_guard+1点恒虚警处理噪声由后面Nref个点决定
    for r_n = 1:N_guard+1
%         if(temp(r_n)<(amd_factor*avg_whole))
%             continue;
%         end
        sum1 = 0;
        for j=1+N_guard+1:N_guard+N_ref
            sum1 = sum1+temp(j);
        end
        
        temp1 = temp(1+N_guard+1:N_guard+N_ref);
        [max_data,max_local] = max(temp1);
        sum1 = sum1-max_data;
        temp1(max_local) = 0;
        
        avg1 = (sum1/(N_ref-1));
        CFAR_1(r_n) = avg1*T;
        
        sum1 = 0;
    end
    
    %第Nguard+1点到第Nref+Nguard-1点恒虚警处理时噪声均值由前后各ref个点共同决定
    for r_n = N_guard+2:N_guard+N_ref
        sum1 = 0;
        for j=1:r_n-N_guard-1
            sum1 = sum1+temp(j);
        end
        for j=r_n+N_guard+1:r_n+N_guard+N_ref
            sum1 = sum1+temp(j);
        end
        
        temp_1(1:r_n-N_guard-1) = temp(1:r_n-N_guard-1);
        temp_1(r_n-N_guard:N_ref+r_n-N_guard-1) = temp(r_n+N_guard+1:r_n+N_guard+N_ref);
        [max_data,max_local] = max(temp_1);
        sum1 = sum1-max_data;
        temp_1(max_local) = 0;
        
        avg1 = (sum1/(N_ref+r_n-N_guard-1-1));
        CFAR_1(r_n) = avg1*T;
        sum1 = 0;
    end
    sum2 = 0;
    
    %正常点的处理
    for r_n=N_ref+N_guard+1:L_v-N_ref-N_guard
        if i<14
            
            for j=r_n-(N_ref+N_guard):r_n-N_guard-1
                sum1 = sum1+temp(j);
            end
            
            for j=r_n+N_guard+1:r_n+(N_ref+N_guard)
                sum1 = sum1+temp(j);
            end
            
            temp_2(1:N_ref) = temp(r_n-(N_ref+N_guard):r_n-N_guard-1);
            temp_2(N_ref+1:2*N_ref) = temp(r_n+N_guard+1:r_n+(N_ref+N_guard));
            [MAX,max_local] = max(temp_2);
            sum1 = sum1-MAX;
            temp_2(max_local) = 0;
            
            [MAX,max_local] = max(temp_2);
            sum1 = sum1-MAX;
            temp_2(max_local) = 0;
            
            [MAX,max_local] = max(temp_2);
            sum1 = sum1-MAX;
            temp_2(max_local) = 0;
            
            [MAX,max_local] = max(temp_2);
            sum1 = sum1-MAX;
            temp_2(max_local) = 0;
            
            [MAX,max_local] = max(temp_2);
            sum1 = sum1-MAX;
            temp_2(max_local) = 0;
            
            avg1 = ((sum1)/(N_ref*2-5));
                
            CFAR_1(r_n) = avg1*T;
            sum1 = 0;
            sum2 = 0;
        end
        
        if i>=14
            for j=r_n-(N_ref+N_guard):r_n-N_guard-1
                sum1 = sum1+temp(j);
            end
            [MAX,max_local] = max(temp(r_n-(N_ref+N_guard):r_n-N_guard-1));
            sum1 = sum1-MAX;
            
            avg1 = ((sum1)/(N_ref-1));
            
            for j=r_n+N_guard+1:r_n+(N_ref+N_guard)
                sum2 = sum2+temp(j);
            end
            
            [MAX,max_local] = max(temp(r_n+N_guard+1:r_n+(N_guard+N_ref)));
            sum2 = sum2-MAX; 
            avg2 = ((sum2)/(N_ref-1));
            
            if(avg1 >= avg2)
                avg1 = avg2;
            else
                avg1 = avg1;
            end
            
            factor_temp = factor;
            
            if ((i>25&&i<50) && (avg1<(avg_whole*0.6)) )
                factor_temp = factor+6;
            end
            if ((i<14) && (avg1>avg_whole*0.9))
                factor_temp  = factor -2;
            end
            
            CFAR_1(r_n) = avg1*T;
            sum1 = 0;
            sum2 = 0;
        end
    end
    
    %第Ns-Nguard-Nref点到第Ns-1-Nguard-1点的噪声由前后Nref个点共同决定
    for r_n = L_v-N_guard-N_ref+1:L_v-N_guard-1
        for j=r_n-(N_ref+N_guard):r_n-N_guard-1
            sum1 = sum1+temp(j);
        end       
        for j = r_n+N_guard+1:L_v
            sum1 = sum1+temp(j);
        end    
        avg1 = (sum1/(N_ref+L_v-r_n-N_guard));
        CFAR_1(r_n) = avg1*T;
        sum1 = 0;
    end
       
    %第Ns-Nguard-1点到点Ns-1点噪声均值由其前面的Nref个点决定
    for r_n = L_v-N_guard:L_v
        for j=r_n-(N_ref+N_guard):r_n-N_guard-1
            sum1 = sum1+temp(j);
        end
        avg1 = (sum1/N_ref);
        CFAR_1(r_n) = avg1*T;
        sum1 = 0;
    end
    
    
    %***************************************************************************************************
    temp=AmbData(i,:);
    avg_whole = mean(temp);
    
    if i<6
        factor = 14;
        amd_factor = 0.5;
    end
    
    if i<35 && i>=6
        factor = 7.5;
        amd_factor = 0.5;
    end
    
    if i<50 && i>=35
        factor = 8.5;
        amd_factor = 0.8;
    end
    
    if i<80 && i>=50
        factor = 8.5;
        amd_factor = 1.2;
    end
    
    if i>=80
        factor = 11;
        amd_factor = 3;
    end
    
    %第1点到第N_guard+1点恒虚警处理噪声由后面Nref个点决定
    for r_n = 1:N_guard+1
%         if(temp(r_n)<(amd_factor*avg_whole))
%             continue;
%         end
        sum1 = 0;
        for j=1+N_guard+1:N_guard+N_ref
            sum1 = sum1+temp(j);
        end
        
        temp1 = temp(1+N_guard+1:N_guard+N_ref);
        [max_data,max_local] = max(temp1);
        sum1 = sum1-max_data;
        temp1(max_local) = 0;
        
        avg1 = (sum1/(N_ref-1));
        CFAR_2(r_n) = avg1*factor;
        sum1 = 0;
    end
    
    %第Nguard+1点到第Nref+Nguard-1点恒虚警处理时噪声均值由前后各ref个点共同决定
    for r_n = N_guard+2:N_guard+N_ref
        sum1 = 0;
        for j=1:r_n-N_guard-1
            sum1 = sum1+temp(j);
        end
        for j=r_n+N_guard+1:r_n+N_guard+N_ref
            sum1 = sum1+temp(j);
        end
        
        temp_1(1:r_n-N_guard-1) = temp(1:r_n-N_guard-1);
        temp_1(r_n-N_guard:N_ref+r_n-N_guard-1) = temp(r_n+N_guard+1:r_n+N_guard+N_ref);
        [max_data,max_local] = max(temp_1);
        sum1 = sum1-max_data;
        temp_1(max_local) = 0;
        
        avg1 = (sum1/(N_ref+r_n-N_guard-1-1));
        CFAR_2(r_n) = avg1*factor;
        sum1 = 0;
    end
    sum2 = 0;
    
    %正常点的处理
    for r_n=N_ref+N_guard+1:L_v-N_ref-N_guard
        if i<14          
            for j=r_n-(N_ref+N_guard):r_n-N_guard-1
                sum1 = sum1+temp(j);
            end
            
            for j=r_n+N_guard+1:r_n+(N_ref+N_guard)
                sum1 = sum1+temp(j);
            end
            
            temp_2(1:N_ref) = temp(r_n-(N_ref+N_guard):r_n-N_guard-1);
            temp_2(N_ref+1:2*N_ref) = temp(r_n+N_guard+1:r_n+(N_ref+N_guard));
            [MAX,max_local] = max(temp_2);
            sum1 = sum1-MAX;
            temp_2(max_local) = 0;
            
            [MAX,max_local] = max(temp_2);
            sum1 = sum1-MAX;
            temp_2(max_local) = 0;
            
            [MAX,max_local] = max(temp_2);
            sum1 = sum1-MAX;
            temp_2(max_local) = 0;
            
            [MAX,max_local] = max(temp_2);
            sum1 = sum1-MAX;
            temp_2(max_local) = 0;
            
            [MAX,max_local] = max(temp_2);
            sum1 = sum1-MAX;
            temp_2(max_local) = 0;
            
            avg1 = ((sum1)/(N_ref*2-5));
            
            factor_temp = factor;
            if ( (i>25&&i<50) && (avg1<(avg_whole*0.6)) )
                factor_temp = factor+6;
            end
            
            if ( (i<14) && (avg1>(avg_whole*0.9)) )
                factor_temp = factor-2;
            end
                CFAR_2(r_n) = avg1*factor;
            sum1 = 0;
            sum2 = 0;
        end
        
        if i>=14           
            for j=r_n-(N_ref+N_guard):r_n-N_guard-1
                sum1 = sum1+temp(j);
            end
            [MAX,max_local] = max(temp(r_n-(N_ref+N_guard):r_n-N_guard-1));
            sum1 = sum1-MAX;
            
            avg1 = ((sum1)/(N_ref-1));
            
            for j=r_n+N_guard+1:r_n+(N_ref+N_guard)
                sum2 = sum2+temp(j);
            end
            
            [MAX,max_local] = max(temp(r_n+N_guard+1:r_n+(N_guard+N_ref)));
            sum2 = sum2-MAX; 
            avg2 = ((sum2)/(N_ref-1));
            
            if(avg1 >= avg2)
                avg1 = avg2;
            else
                avg1 = avg1;
            end
            
            factor_temp = factor;
            
            if ((i>25&&i<50) && (avg1<(avg_whole*0.6)) )
                factor_temp = factor+6;
            end
            if ((i<14) && (avg1>avg_whole*0.9))
                factor_temp  = factor -2;
            end
            CFAR_2(r_n) = avg1*factor_temp;
            sum1 = 0;
            sum2 = 0;
        end
    end
    
    %第Ns-Nguard-Nref点到第Ns-1-Nguard-1点的噪声由前后Nref个点共同决定
    for r_n = L_v-N_guard-N_ref+1:L_v-N_guard-1
        for j=r_n-(N_ref+N_guard):r_n-N_guard-1
            sum1 = sum1+temp(j);
        end       
        for j = r_n+N_guard+1:L_v
            sum1 = sum1+temp(j);
        end    
        avg1 = (sum1/(N_ref+L_v-r_n-N_guard));
        CFAR_2(r_n) = avg1*factor;
        sum1 = 0;
    end
       
    %第Ns-Nguard-1点到点Ns-1点噪声均值由其前面的Nref个点决定
    for r_n = L_v-N_guard:L_v
        for j=r_n-(N_ref+N_guard):r_n-N_guard-1
            sum1 = sum1+temp(j);
        end
        avg1 = (sum1/N_ref);
        CFAR_2(r_n) = avg1*factor;
        sum1 = 0;
    end 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OS_CFAR算法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp=AmbData(i,:);
    OS_CFAR_Factor = 7;
    %第1点到第N_guard+1点恒虚警处理噪声由后面Nref个点决定
    for r_n = 1:N_guard+1
        
        os_index = 1;
        for j=1+N_guard+1:N_guard+N_ref*2
            temp_OS_1(os_index) = temp(j);
            os_index = os_index + 1;
        end
        
        for kk=1:os_index-1
            [max_os,index_os] = max(temp_OS_1);
            temp_OS_1(index_os) = 0;
            temp_OS_2(os_index-kk) = max_os;
        end
        
        temp_OS_2(os_index-1) = 0;
        temp_OS_2(os_index-2) = 0;
        os_mean = sum(temp_OS_2)/(os_index-3);
        
        OS_CFAR_gate(r_n) = os_mean*OS_CFAR_Factor;   
    end
    
    %第Nguard+1点到第Nref+Nguard-1点恒虚警处理时噪声均值由前后各ref个点共同决定
    for r_n = N_guard+2:N_guard+N_ref
        os_index = 1;
        for j=1:r_n-N_guard-1
            temp_OS_1(os_index) = temp(j);
            os_index = os_index + 1;
        end
        for j=r_n+N_guard+1:r_n+N_guard+N_ref
            temp_OS_1(os_index) = temp(j);
            os_index = os_index + 1;
        end
        
        for kk=1:os_index-1
            [max_os,index_os] = max(temp_OS_1);
            temp_OS_1(index_os) = 0;
            temp_OS_2(os_index-kk) = max_os;
        end
        
        temp_OS_2(os_index-1) = 0;
        temp_OS_2(os_index-2) = 0;
        os_mean = sum(temp_OS_2)/(os_index-3);
        
        OS_CFAR_gate(r_n) = os_mean*OS_CFAR_Factor; 
    end
    
    %正常点的处理
    for r_n=N_ref+N_guard+1:L_v-N_ref-N_guard
        os_index = 1;   
        for j=r_n-(N_ref+N_guard):r_n-N_guard-1
            temp_OS_1(os_index) = temp(j);
            os_index = os_index + 1;
        end

        for j=r_n+N_guard+1:r_n+(N_ref+N_guard)
            temp_OS_1(os_index) = temp(j);
            os_index = os_index + 1;
        end

        for kk=1:os_index-1
            [max_os,index_os] = max(temp_OS_1);
            temp_OS_1(index_os) = 0;
            temp_OS_2(os_index-kk) = max_os;
        end
        
        temp_OS_2(os_index-1) = 0;
        temp_OS_2(os_index-2) = 0;
        os_mean = sum(temp_OS_2)/(os_index-3);
        
        OS_CFAR_gate(r_n) = os_mean*OS_CFAR_Factor; 
    end
    
    %第Ns-Nguard-Nref点到第Ns-1-Nguard-1点的噪声由前后Nref个点共同决定
    for r_n = L_v-N_guard-N_ref+1:L_v-N_guard-1
        os_index = 1;
        for j=r_n-(N_ref+N_guard):r_n-N_guard-1
            temp_OS_1(os_index) = temp(j);
            os_index = os_index + 1;
        end       
        for j = r_n+N_guard+1:L_v
            temp_OS_1(os_index) = temp(j);
            os_index = os_index + 1;
        end    
        for kk=1:os_index-1
            [max_os,index_os] = max(temp_OS_1);
            temp_OS_1(index_os) = 0;
            temp_OS_2(os_index-kk) = max_os;
        end
        temp_OS_2(os_index-1) = 0;
        temp_OS_2(os_index-2) = 0;
        os_mean = sum(temp_OS_2)/(os_index-3);
        
        OS_CFAR_gate(r_n) = os_mean*OS_CFAR_Factor; 
    end
       
    %第Ns-Nguard-1点到点Ns-1点噪声均值由其前面的Nref个点决定
    for r_n = L_v-N_guard:L_v
        os_index = 1;
        for j=r_n-(N_ref*2+N_guard):r_n-N_guard-1
            temp_OS_1(os_index) = temp(j);
            os_index = os_index + 1;
        end
        
        for kk=1:os_index-1
            [max_os,index_os] = max(temp_OS_1);
            temp_OS_1(index_os) = 0;
            temp_OS_2(os_index-kk) = max_os;
        end
        temp_OS_2(os_index-1) = 0;
        temp_OS_2(os_index-2) = 0;
        os_mean = sum(temp_OS_2)/(os_index-3);
        
        OS_CFAR_gate(r_n) = os_mean*OS_CFAR_Factor; 
    end 
    
    
    plot(temp(1:64),'-b*');
    hold on;
    plot(CFAR_1(1:64),'-.m');
    hold on;
    plot(CFAR_2(1:64),'-.r');
    hold on;
    plot(OS_CFAR_gate(1:64),'-.black');
    hold on;
    legend('Raw Data','CFAR 1','CFAR 2','OS CFAR');
    pause(0.3);
    clf
    
    
end