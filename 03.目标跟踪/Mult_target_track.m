
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����24G�״�CANԭʼ��Ϣ�Ķ�Ŀ����ٷ����㷨
%20180601��ʵ�ֻ���������
%20180609��1��������򷨹��������͵㼣
%          2��������Ŀ������˲���Ϊ�ɿ������˲���ά��           
%          3���ɿ���������ʼ����ֱ�۷�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear all;
close all;

T = 0.05;

global Pd Pg gamma Q R ffun hfun;

Pd=1;       %������ 
Pg=0.99;    %��ȷ��������������ڵø��� 
gamma = 1; 
    
% ��������Э����
Q =[1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];

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


R11=1; R22=1; R12=0; R21=0;

%��ʼЭ���� 
P=[R11 R11/T R12 R12/T; R11/T 2*R11/T^2 R12/T 2*R12/T^2; 
   R21 R21/T R22 R22/T; R21/T 2*R21/T^2 R22/T 2*R22/T^2];   

global confirm_track_ID temp_track_ID confirm_target_number temp_target_number Conv_P middle_track_ID middle_target_number;
Target_Table1 = struct( 'status',{},...
                        'counter',{},...
                        'no_related',{},...
                        'M',zeros(4,1),...
                        'X',zeros(4,1),...
                        'Z',zeros(2,1),...
                        'P',zeros(4,4),...
                        'FROM',{},...
                        'Middle_counter',0,...
                        'kekao_counter',0,...
                        'X_AVE',{} );
Conv_P = P;
for i=1:256
    Target_Table1(i).status = 0;
    Target_Table1(i).counter = 0;
    Target_Table1(i).no_related = 0;
    Target_Table1(i).M = zeros(4,1);
    Target_Table1(i).X = zeros(4,1);
    Target_Table1(i).Z = zeros(2,1);
    Target_Table1(i).P = zeros(4,4);
    Target_Table1(i).FROM = 0;
    Target_Table1(i).X_AVE = 0;
end

global Track_Table;
Track_Table = Target_Table1;
                       
confirm_track_ID = zeros(64,1);
temp_track_ID = zeros(64,1);
middle_track_ID = zeros(64,1);
confirm_target_number = 0;
temp_target_number = 0;
middle_target_number = 0;

fid = fopen('Target_info.txt','r');
data = textscan(fid,'%f %f %f');



line_index = 1;
target_info = zeros(64,4);
time = 0;
while 1
    clf
    time = time + 1;
    target_num = data{1}(line_index,1);
    if target_num>0 
        for j = 1:target_num
            
            RR = data{1}(line_index+j,1); %R
            V  = data{2}(line_index+j,1); %V
            A  = data{3}(line_index+j,1)+25; %A

            target_info(j,1) = RR*sin(abs(A)*3.14/180);  %X
            target_info(j,3) = RR*cos(abs(A)*3.14/180);  %Y
            target_info(j,2) = 0.01;%V*sin(A*3.14/180);        %VX
            target_info(j,4) = V*cos(A*3.14/180);        %VY
           
            if A < 0
                target_info(j,1) = -1*target_info(j,1);
            end
        end
    end

    % �㼣Ԥ����
%     fprintf('*********************** step %d ***********************',time)
    % ���ԭʼĿ����Ϣ
    ori_target_info = target_info(1:target_num,:);
    [target_info1,target_num1] = Pre_Deal(target_info,target_num);
   
    % ���Ԥ�����ԭʼĿ����Ϣ
    before_info = target_info1(1:target_num1,:);
    x = target_info1(1:target_num1,1);y = target_info1(1:target_num1,3); plot(x,y,'k*');hold on   
    %����������
    [new_target_number,new_target_info, new_target_number_1,new_target_info_1,] = track_func(target_num1,target_info1);
    % ������ɿ��������������˲���Ŀ����Ϣ
    kalman_info = new_target_info(1:new_target_number,:);
    x = new_target_info(1:new_target_number,1);    y = new_target_info(1:new_target_number,3); plot(x,y,'r.','MarkerSize',15);    hold on 
    % �����ʱ����Ŀ����Ϣ
    linear_info = new_target_info_1(1:new_target_number_1,:);
    x = new_target_info_1(1:new_target_number_1,1);    y = new_target_info_1(1:new_target_number_1,3); plot(x,y,'g.','MarkerSize',15);    hold on
    axis([-40 40 0 70]); grid on; 
    line_index = line_index + target_num +1;
    pause(0.05);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �״�����Ԥ������
% ���ٽ��������㼣����ͬһ��Ŀ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out,num_out] = Pre_Deal(info,num)
    for i = 1:num-1
        delta_x = abs(info(i,1) - info(i+1,1));
        delta_y = abs(info(i,3) - info(i+1,3));

        if delta_x<1 && delta_y<1
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