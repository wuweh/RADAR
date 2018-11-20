%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   航迹管理函数说明：
%   input:  1）N_target：             目标个数
%           2）Target_Info：          目标信息
%   output: 1）Target_Number：        可靠航迹个数
%           2）Target_Info_Output：   可靠航迹信息
%           3）Target_Number_1：      临时航迹个数
%           4）Target_Info_Output_1： 临时航迹信息
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%   Change_History:
%       1）20180609：利用卡方分布设置门限，对落入门限内的量测点采用最近领域法找关联点，滤波采用4维参数标准卡尔曼滤波
%       2）20180610：修改可靠航迹点迹相关方法，直接用XY轴方向门限判定点迹是否与航迹相关
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [Target_Number, Target_Info_Output,Target_Number_1, Target_Info_Output_1] = track_func(N_target,Target_Info)
    global Track_Table confirm_track_ID temp_track_ID confirm_target_number temp_target_number Conv_P middle_track_ID middle_target_number;    
      
    % 存在确认航迹
    if confirm_target_number>0 
        for i=1:confirm_target_number
            deleta_X = abs(Target_Info(1:N_target,1) - Track_Table(confirm_track_ID(i)).X(1));
            deleta_Y = abs(Target_Info(1:N_target,3) - Track_Table(confirm_track_ID(i)).X(3));
            d_2 = (deleta_X.^2 + deleta_Y.^2)/4;
            [value, index] = min(d_2);

            if value<3
                n = 1;
            else
                n = 0;
            end

            if n>0
                 % 删除已关联的点迹信息
                Track_Table(confirm_track_ID(i)).M = Target_Info(index,:)';
                Track_Table(confirm_track_ID(i)).FROM = Track_Table(confirm_track_ID(i)).FROM + 1;
                
                %删除已关联的点迹信息
                for k=index:N_target-1 
                    Target_Info(k,:) = Target_Info(k+1,:);
                end         
                N_target = N_target - 1; 
            
                %开始卡尔曼滤波
                x = Track_Table(confirm_track_ID(i)).X;
                y = Track_Table(confirm_track_ID(i)).M;
                p = Track_Table(confirm_track_ID(i)).P;               
                [x1,p1] = kalman(x,y,p);                
                Track_Table(confirm_track_ID(i)).X = x1;
                Track_Table(confirm_track_ID(i)).P = p1;
            else 
                 if  Track_Table(confirm_track_ID(i)).no_related == 4 %连续四次未关联上目标
                    %清空本航迹所有信息
                    Track_Table(confirm_track_ID(i)).status = 0;%delete
                    Track_Table(confirm_track_ID(i)).M = zeros(4,1);
                    Track_Table(confirm_track_ID(i)).X = zeros(4,1);
                    Track_Table(confirm_track_ID(i)).Z = zeros(2,1);
                    Track_Table(confirm_track_ID(i)).counter = 0;
                    Track_Table(confirm_track_ID(i)).no_related = 0;
                    Track_Table(confirm_track_ID(i)).Conv = zeros(6,6);
                    Track_Table(confirm_track_ID(i)).FROM = 0;
                    Track_Table(confirm_track_ID(i)).X_AVE = 0;
                else
                    % 线性预测
                    Track_Table(confirm_track_ID(i)).no_related = Track_Table(confirm_track_ID(i)).no_related + 1;
                    Track_Table(confirm_track_ID(i)).M = zeros(4,1);
                    Track_Table(confirm_track_ID(i)).Z(1) = 0;
                    Track_Table(confirm_track_ID(i)).Z(2) = 0; %本次未关联上目标，用上周期预测值作为本次量测值进行线性预测
                    [Track_Table(confirm_track_ID(i)).X] = linear_filter(Track_Table(confirm_track_ID(i)).X); 
                    Track_Table(confirm_track_ID(i)).FROM = Track_Table(confirm_track_ID(i)).FROM + 1;
                end
            end          
       end
    end
    
       % 存在中间航迹
     if middle_target_number>0 
        for i=1:middle_target_number
                n = 0;
                deleta_X = abs(Target_Info(1:N_target,1) - Track_Table(middle_track_ID(i)).X(1));
                deleta_Y = abs(Target_Info(1:N_target,3) - Track_Table(middle_track_ID(i)).X(3));
                d_2 = (deleta_X.^2 + deleta_Y.^2)/4;
                [value, index] = min(d_2);
                if value<7
                     n = 1;
                     
                     if value < 3
                       Track_Table(middle_track_ID(i)).kekao_counter = Track_Table(middle_track_ID(i)).kekao_counter + 1;
                     else
                        if Track_Table(middle_track_ID(i)).kekao_counter ~= 0
                            Track_Table(middle_track_ID(i)).kekao_counter =Track_Table(middle_track_ID(i)).kekao_counter - 1;
                        end
                     end
                     
                else
                    n = 0;
                end       
                
            if n==1
                Track_Table(middle_track_ID(i)).no_related = 0;
                Track_Table(middle_track_ID(i)).M = Target_Info(index,:)';
                Track_Table(middle_track_ID(i)).counter = Track_Table(middle_track_ID(i)).counter + 1;
                  for k=index:N_target-1 %删除已关联的点迹信息
                      Target_Info(k,:) = Target_Info(k+1,:);
                  end         
                  N_target = N_target - 1; 
                  
                  if Track_Table(middle_track_ID(i)).kekao_counter == 2  % 连续两次，升级为可靠航迹
                            Track_Table(middle_track_ID(i)).status = 3;  
                            Track_Table(middle_track_ID(i)).kekao_counter  = 0;
                            Track_Table(middle_track_ID(i)).Middle_counter = 0;
                            Track_Table(middle_track_ID(i)).no_related = 0;
                            
                            Track_Table(middle_track_ID(i)).X = Track_Table(middle_track_ID(i)).M;
                            Track_Table(middle_track_ID(i)).P= Conv_P;
                            Track_Table(middle_track_ID(i)).FROM = 1;
                            %开始卡尔曼滤波
                            x = Track_Table(middle_track_ID(i)).X;
                            y = Track_Table(middle_track_ID(i)).M;
                            p = Track_Table(middle_track_ID(i)).P;
                            [Track_Table(middle_track_ID(i)).X, Track_Table(middle_track_ID(i)).P] = kalman(x,y,p);

                  else
                            [Track_Table(middle_track_ID(i)).X] = linear_filter(Track_Table(middle_track_ID(i)).X);
                   end        
            else
                Track_Table(middle_track_ID(i)).no_related  = Track_Table(middle_track_ID(i)).no_related + 1;
                Track_Table(middle_track_ID(i)).counter = 0;
                if Track_Table(middle_track_ID(i)).no_related == 2 %连续两次没有关联点迹，删除中间航迹
                    Track_Table(middle_track_ID(i)).status = 0;%delete
                    Track_Table(middle_track_ID(i)).M = zeros(4,1);
                    Track_Table(middle_track_ID(i)).X = zeros(4,1);
                    Track_Table(middle_track_ID(i)).Z = zeros(2,1);
                    Track_Table(middle_track_ID(i)).Conv = zeros(6,6);
                    Track_Table(middle_track_ID(i)).counter = 0;
                    Track_Table(middle_track_ID(i)).no_related = 0;
                    Track_Table(middle_track_ID(i)).FROM = 0;
                    Track_Table(middle_track_ID(i)).kekao_counter  = 0;
                    Track_Table(middle_track_ID(i)).Middle_counter = 0;
                else
                    if Track_Table(middle_track_ID(i)).kekao_counter ~= 0
                        Track_Table(middle_track_ID(i)).kekao_counter = Track_Table(middle_track_ID(i)).kekao_counter - 1;
                    end
                        % 线性预测
                        [Track_Table(middle_track_ID(i)).X] = linear_filter(Track_Table(middle_track_ID(i)).X);
                end
            end  
        end
     end
     
    % 存在临时航迹
     if temp_target_number>0 %存在临时航迹
        for i=1:temp_target_number
                deleta_X = abs(Target_Info(1:N_target,1) - Track_Table(temp_track_ID(i)).X(1));
                deleta_Y = abs(Target_Info(1:N_target,3) - Track_Table(temp_track_ID(i)).X(3));
                d_2 = (deleta_X.^2 + deleta_Y.^2)/4;
                [value, index] = min(d_2);
                if value<10
                    n = 1;
                    if value < 3
                       Track_Table(temp_track_ID(i)).kekao_counter = Track_Table(temp_track_ID(i)).kekao_counter + 1;
                    else
                        Track_Table(temp_track_ID(i)).Middle_counter = Track_Table(temp_track_ID(i)).Middle_counter + 1;
                    end
                else
                    n = 0;
                end
                   
            if n==1
                  Track_Table(temp_track_ID(i)).M = Target_Info(index,:)';
                  Track_Table(temp_track_ID(i)).counter = Track_Table(temp_track_ID(i)).counter + 1;
                  for k=index:N_target-1 %删除已关联的点迹信息
                      Target_Info(k,:) = Target_Info(k+1,:);
                  end         
                  N_target = N_target - 1; 
                  if Track_Table(temp_track_ID(i)).counter == 5
                      %此处判断升级为可靠航迹还是中间航迹
                      if Track_Table(temp_track_ID(i)).kekao_counter >= 3
                            Track_Table(temp_track_ID(i)).status = 3;  %升级为可靠航迹
                            Track_Table(temp_track_ID(i)).X = Track_Table(temp_track_ID(i)).M;
                            Track_Table(temp_track_ID(i)).P= Conv_P;
                            Track_Table(temp_track_ID(i)).FROM = 1;
                            %开始卡尔曼滤波
                            x = Track_Table(temp_track_ID(i)).X;
                            y = Track_Table(temp_track_ID(i)).M;
                            p = Track_Table(temp_track_ID(i)).P;
                            [x1,p1] = kalman(x,y,p);
                            Track_Table(temp_track_ID(i)).X = x1;
                            Track_Table(temp_track_ID(i)).P = p1; 
                            Track_Table(temp_track_ID(i)).kekao_counter  = 0;
                            Track_Table(temp_track_ID(i)).Middle_counter = 0;
                      else
                            Track_Table(temp_track_ID(i)).status = 2;   %升级为中间航迹
                            y = Track_Table(temp_track_ID(i)).M;
                            [x1] = linear_filter(y);
                            Track_Table(temp_track_ID(i)).X = x1;
                            Track_Table(temp_track_ID(i)).kekao_counter  = 0;
                            Track_Table(temp_track_ID(i)).Middle_counter = 0;
                      end        
                  else
                        % 线性预测
                    [Track_Table(temp_track_ID(i)).X] = linear_filter(Track_Table(temp_track_ID(i)).X);
                  end
            end
            if n == 0
                Track_Table(temp_track_ID(i)).status = 0;%delete
                Track_Table(temp_track_ID(i)).M = zeros(4,1);
                Track_Table(temp_track_ID(i)).X = zeros(4,1);
                Track_Table(temp_track_ID(i)).Z = zeros(2,1);
                Track_Table(temp_track_ID(i)).Conv = zeros(6,6);
                Track_Table(temp_track_ID(i)).counter = 0;
                Track_Table(temp_track_ID(i)).no_related = 0;
                Track_Table(temp_track_ID(i)).FROM = 0;
            end  
        end
     end

    % 没有临时航迹也没有确认航迹
    if confirm_target_number == 0
        if temp_target_number==0
            if N_target>0
                for i=1:N_target
                    Track_Table(i).status = 1; %1表示临时航迹
                    Track_Table(i).M = Target_Info(i,:)';
                    Track_Table(i).counter = 1;
                    temp_target_number = temp_target_number + 1;
                    % 线性预测
                    [Track_Table(i).X] = linear_filter(Track_Table(i).M); 
                end
                N_target = 0;
            end
        end
    end
    
    % 本次量测未关联上的点
    if N_target>0
        for i=1:N_target   
            for j=1:64
                if Track_Table(j).status == 0 %该位置可用
                    Track_Table(j).status = 1;
                    Track_Table(j).M = Target_Info(i,:)';
                    Track_Table(j).counter = 1;
                    % 线性预测
                    [Track_Table(j).X] = linear_filter(Track_Table(j).M); 
                    break;
                end
            end
        end
    end
    N_target = 0;
    
    % 整理航迹序号
    temp_target_number = 0; 
    temp_track_ID = zeros(256,1);
    middle_target_number = 0;
    middle_track_ID = zeros(256,1);
    confirm_target_number = 0;
    confirm_track_ID = zeros(256,1);
    
    for ii=1:256
        if  Track_Table(ii).status == 1
            temp_target_number = temp_target_number+ 1;
            temp_track_ID(temp_target_number) = ii;       
        end

        if  Track_Table(ii).status == 2
            middle_target_number = middle_target_number+ 1;
            middle_track_ID(middle_target_number) = ii;       
        end

        if  Track_Table(ii).status == 3
            confirm_target_number = confirm_target_number+ 1;
            confirm_track_ID(confirm_target_number) = ii;       
        end
    end
    
    % 输出结果
    Target_Number = confirm_target_number;
    Target_Info_Output = zeros(Target_Number,4);
    for i=1:Target_Number
        Target_Info_Output(i,1:4) = Track_Table(confirm_track_ID(i)).X';
        %对X信息进行平滑
        counters = Track_Table(confirm_track_ID(i)).FROM;
        x_value = Track_Table(confirm_track_ID(i)).X_AVE;
        Target_Info_Output(i,1) = (((counters -1)*x_value) + Target_Info_Output(i,1))/counters;
        Track_Table(confirm_track_ID(i)).X_AVE = Target_Info_Output(i,1);
        Target_Info_Output(i,5:8) = Track_Table(confirm_track_ID(i)).M;
    end
    
    Target_Number_1 = temp_target_number;
    Target_Info_Output_1 = zeros(Target_Number_1,4);
    for i=1:Target_Number_1
        Target_Info_Output_1(i,1:4) = Track_Table(temp_track_ID(i)).X';
        Target_Info_Output_1(i,5) = Track_Table(temp_track_ID(i)).counter;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   卡尔曼滤波函数说明：
%   input:  1）X：    上一次预测值
%           2）Z：    量测值
%           3）P：    协方差矩阵
%   output: 1）x：    当前预测值
%           2）p：    后验协方差矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [x,p] = kalman(X,Z,P)
    global Q R ffun hfun;     
    n=numel(X);
    x1 = ffun*X;
    P1 = ffun*P*ffun'+Q;
    z = hfun*x1;
    S = hfun*P1*hfun'+R;
    gain = P1*hfun'/S;
    
    x = x1+gain*(Z-z);
    p = (eye(n)-gain*hfun)*P1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   线性滤波函数说明：
%   input:  1）x：当前状态
%   output: 1）X：预测状态
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X] = linear_filter(x)  
    global ffun;
    X = ffun*x;
end

            
