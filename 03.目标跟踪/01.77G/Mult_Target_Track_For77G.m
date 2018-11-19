
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%基于24G雷达CAN原始信息的多目标跟踪仿真算法
%20180601：实现基本框架设计
%20180609：1）最近邻域法关联航迹和点迹
%          2）最基本的卡尔曼滤波作为可靠航迹滤波四维）           
%          3）可靠航迹的起始采用直观法
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear all;
close all;

load('INFO_DATA.mat');

T = 0.1;

global Pd Pg gamma Q R ffun hfun;

Pd=1;       %检测概率 
Pg=0.99;    %正确量测落入跟踪门内得概率 
gamma = 1; 
    
% 过程噪声协方差
Q =[0.25 0 0 0;
    0 0.25 0 0;
    0 0 0.25 0;
    0 0 0 0.25];

R =[1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];

ffun = [1 T 0 0;
        0 1 0 0;
        0 0 1 T;
        0 0 0 1];

hfun = [1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        0 0 0 1];


R11=0.01; R22=0.01; R12=0; R21=0;

%初始协方差 
P=[R11 R11/T R12 R12/T; R11/T 2*R11/T^2 R12/T 2*R12/T^2; 
   R21 R21/T R22 R22/T; R21/T 2*R21/T^2 R22/T 2*R22/T^2];   

global confirm_track_ID temp_track_ID confirm_target_number temp_target_number Conv_P;
Target_Table1 = struct( 'status',{},...
                        'counter',{},...
                        'no_related',{},...
                        'M',zeros(4,1),...
                        'X',zeros(4,1),...
                         'X_out',zeros(4,1),...
                        'Z',zeros(2,1),...
                        'P',zeros(4,4),...
                        'FROM',{},...
                        'X_AVE',{} );
Conv_P = P;
for i=1:256
    Target_Table1(i).status = 0;
    Target_Table1(i).counter = 0;
    Target_Table1(i).no_related = 0;
    Target_Table1(i).M = zeros(4,1);
    Target_Table1(i).X = zeros(4,1);
    Target_Table1(i).X_out = zeros(4,1);
    Target_Table1(i).Z = zeros(2,1);
    Target_Table1(i).P = zeros(4,4);
    Target_Table1(i).FROM = 0;
    Target_Table1(i).X_AVE = 0;
end

global Track_Table;
Track_Table = Target_Table1;
                       
confirm_track_ID = zeros(64,1);
temp_track_ID = zeros(64,1);
confirm_target_number = 0;
temp_target_number = 0;

target_info = zeros(64,4);
time = 0;
target_num = 0;
index = 1;
counter = 0;

while 1   
    counter = counter + 1;
    a = INFO_DATA{1,1}(index,1);   
    if a > 1000000
        time = a+18;
        index = index+1;
        
        time = time-20180621000000000;
        if time == 103640224
            pause(0.1);
        end
        step = sprintf('*********************** step %d ***********************',time) 
        time = num2str(time);
         
        if target_num ~= 0
            % 输出原始目标信息
            [target_info1,target_num1] = Pre_Deal(target_info,target_num);    
            figure(1); clf
            plot(target_info1(1:target_num1,1),target_info1(1:target_num1,3),'b.','MarkerSize',20);hold on   
            X_range = sprintf('%.2f  ', target_info1(1:target_num1,1))
            Y_range = sprintf('%.2f  ', target_info1(1:target_num1,3))
            v_y = sprintf('%.2f  ', target_info1(1:target_num1,4))
        else
            target_num1 =0;
            target_info1 = zeros(1,4);
        end
            %航迹管理函数
            [new_target_number,new_target_info, new_target_number_1,new_target_info_1,] = track_func(target_num1,target_info1);
            % 输出（可靠航迹）卡尔曼滤波后目标信息
%             new_target_info
            plot(new_target_info(1:new_target_number,1),new_target_info(1:new_target_number,3),'r.','MarkerSize',20);hold on   
            plot(new_target_info_1(1:new_target_number_1,1),new_target_info_1(1:new_target_number_1,3),'g.','MarkerSize',20);hold on   
            
            axis([-20 20 0 200]);  grid on     
            title(time);
            pause(0.1);
            target_num = 0;
    else
        %提取目标信息
        target_num = target_num +1;
        Range = INFO_DATA{1,1}(index,1); %R
        
        index=index+1;
        A  = INFO_DATA{1,1}(index,1); %A
        
        index=index+1;
        V  = INFO_DATA{1,1}(index,1); %V     

        target_info(target_num,1) = Range*sin(A*3.14/180);  %X
        target_info(target_num,3) = Range*cos(A*3.14/180);  %Y
        target_info(target_num,2) = 0.01;        %VX
        target_info(target_num,4) = V;        %VY      
        
        if A < 0
           target_info(target_num,1) = -1*target_info(target_num,1);
        end
        index = index+4;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 雷达数据预处理函数
% 将临近的两个点迹视作同一个目标
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out,num_out] = Pre_Deal(info,num)
    if num == 1
        num_out = 1;
        out = info;
    else
        for i = 1:num-1
            delta_x = abs(info(i,1) - info(i+1,1));
            delta_y = abs(info(i,3) - info(i+1,3));

            if delta_x<3 && delta_y<2
                info(i+1,1) = (info(i,1) + info(i+1,1))/2;
                info(i+1,3) = (info(i,3) + info(i+1,3))/2;
                info(i,:) = zeros(1,4);
            end   
        end

        num_out=0;
        for i=1:num
            if info(i,1) ~= 0
                num_out = num_out+1;
                out(num_out,:) = info(i,:);
            end
        end
    end  
end