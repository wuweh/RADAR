%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   航迹管理函数说明：
%   input:  1）N_target：             目标个数
%               2）Target_Info：          目标信息
%   output: 1）Target_Number：        可靠航迹个数
%           2）Target_Info_Output：   可靠航迹信息
%           3）Target_Number_1：      临时航迹个数
%           4）Target_Info_Output_1： 临时航迹信息
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%   Change_History:
%       1）20180609：利用卡方分布设置门限，对落入门限内的量测点采用最近领域法找关联点，滤波采用4维参数标准卡尔曼滤波
%       2）20180610：修改可靠航迹点迹相关方法，直接用XY轴门限判定点迹是否与航迹相关
%       3) 20180710：修正可靠航迹再次关联到点迹时，没有清空未关联次数的bug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [Target_Number, Target_Info_Output,Target_Number_1, Target_Info_Output_1] = track_func(N_target,Target_Info)
    global Track_Table confirm_track_ID temp_track_ID confirm_target_number temp_target_number Conv_P ffun hfun Q R;    
      
    % 存在确认航迹
    if confirm_target_number>0 
        for i=1:confirm_target_number
            n = 0;
            temp_X = Track_Table(confirm_track_ID(i)).X; 
            temp_P = Track_Table(confirm_track_ID(i)).P;
            
            temp_x1 = ffun*temp_X;  %对上一次的后验值进行一步预测
            temp_P1 = ffun*temp_P*ffun'+Q;
            temp_z = hfun*temp_x1;
            temp_S = hfun*temp_P1*hfun'+R;
            temp_K = temp_P1*hfun'*inv(temp_S); %增益 
            
            %计算卡方值
            if N_target>0
                    temp_d(:,1:N_target) = (Target_Info(1:N_target,:)'-temp_z);
                    for j = 1:N_target
                        temp_d_1(j) = temp_d(:,j)'*inv(temp_S)*temp_d(:,j);
                    end

                    %选择卡方值最小的目标
                    [value,index] = min( temp_d_1(1:N_target));
                    if value < 5 %卡方门限
                        n = 1;
                    end
            end
            
            if n>0
                        Track_Table(confirm_track_ID(i)).M = Target_Info(index,:)';
                        Track_Table(confirm_track_ID(i)).FROM = Track_Table(confirm_track_ID(i)).FROM + 1;
                        Track_Table(confirm_track_ID(i)).no_related = 0;

                        %删除已关联的点迹信息
                        if index<N_target
                            Target_Info(index:N_target-1,:) = Target_Info(index+1:N_target,:);     
                        end
                        N_target = N_target - 1;

                        %开始卡尔曼滤波
                        y = Track_Table(confirm_track_ID(i)).M;       
                        x = temp_x1+temp_K*(y-temp_z);
                        p = (eye(4)-temp_K*hfun)*temp_P1;               
                        Track_Table(confirm_track_ID(i)).X = x;
                        Track_Table(confirm_track_ID(i)).P = p;
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
                                    [Track_Table(confirm_track_ID(i)).X] = linear_filter(Track_Table(confirm_track_ID(i)).X); 
                                    Track_Table(confirm_track_ID(i)).FROM = Track_Table(confirm_track_ID(i)).FROM + 1;
                                    Track_Table(confirm_track_ID(i)).P = Conv_P;
                         end                
            end          
       end
    end
    
    % 存在临时航迹
     if temp_target_number>0 
            for i=1:temp_target_number
                 n = 0;   
                for j = 1:N_target
                    deleta_X = abs(Target_Info(j,1) - Track_Table(temp_track_ID(i)).X(1));
                    deleta_Y = abs(Target_Info(j,3) - Track_Table(temp_track_ID(i)).X(3));

                    if deleta_X<5 && deleta_Y<3 %落入门限范围内
                        n = n + 1;
                        temp(n) = deleta_X^2 + deleta_Y^2;
                        temp_index(n) = j;
                    end
                    if n > 0
                        [~,index] = min(temp(1:n));
                    end
                end

                if n>0
                        Track_Table(temp_track_ID(i)).M = Target_Info(index,:)';
                        Track_Table(temp_track_ID(i)).counter = Track_Table(temp_track_ID(i)).counter + 1;

                        %删除已关联的点迹信息
                        if index<N_target
                            Target_Info(index:N_target-1,:) = Target_Info(index+1:N_target,:);     
                        end
                        N_target = N_target - 1;

                          if Track_Table(temp_track_ID(i)).counter == 4 %连续四次关联成功，升级为可靠航迹
                                    Track_Table(temp_track_ID(i)).status = 2;
                                    Track_Table(temp_track_ID(i)).P= Conv_P;
                                    Track_Table(temp_track_ID(i)).FROM= 1;
                                    Track_Table(temp_track_ID(i)).X = Track_Table(temp_track_ID(i)).M;       
                          else
                                    % 线性预测
                                    y = Track_Table(temp_track_ID(i)).M;     
                                    [x1] = linear_filter(y);
                                    Track_Table(temp_track_ID(i)).X = x1;
                          end
                else
                        Track_Table(temp_track_ID(i)).status = 0;%delete
                        Track_Table(temp_track_ID(i)).M = zeros(4,1);
                        Track_Table(temp_track_ID(i)).X = zeros(4,1);
                        Track_Table(temp_track_ID(i)).X_out= zeros(4,1);
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
                    Track_Table(i).counter = 0;
                    Track_Table(i).P = Conv_P;
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
                    Track_Table(j).counter = 0;
                
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
    temp_track_ID = zeros(64,1);
    confirm_target_number = 0;
    confirm_track_ID = zeros(64,1);
    
    for ii=1:64
       if  Track_Table(ii).status == 1
           temp_target_number = temp_target_number+ 1;
            temp_track_ID(temp_target_number) = ii;       
       end
       
       if  Track_Table(ii).status == 2
           confirm_target_number = confirm_target_number+ 1;
            confirm_track_ID(confirm_target_number) = ii;       
       end
    end
    
    % 输出结果
    Target_Number = confirm_target_number;
    Target_Info_Output = zeros(Target_Number,8);
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
    Target_Info_Output_1 = zeros(Target_Number_1,5);
    for i=1:Target_Number_1
        Target_Info_Output_1(i,1:4) = Track_Table(temp_track_ID(i)).X_out';
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

            
