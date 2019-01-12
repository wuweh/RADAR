clear all
close all
clc

 %%%%%%%%%% ?D_S算法设置 ? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%?
num_Sensor=3;  %传感器数目?
num_Object=6;  %识别对象数目 ?四个目标(机载雷达: ?火控 瞄准 导航 预警)+全集+空集?
num_Period=3;  %测量周期数?

Info=zeros(num_Sensor,num_Object,num_Period);  % 一次判决所需要的信息?
Info(:,:,1)=[0.30 0.40 0.15 0.00 0.15 0.00;
0.30 0.50 0.10 0.00 0.10 0.00;
0.30 0.30 0.20 0.00 0.20 0.00;];

Info(:,:,2)=[0.40 0.20 0.20 0.00 0.20 0.00;
0.50 0.20 0.20 0.00 0.10 0.00;
0.50 0.30 0.10 0.00 0.10 0.00;];

 Info(:,:,3)=[0.50 0.20 0.15 0.00 0.15 0.00;
0.40 0.30 0.10 0.00 0.20 0.00;
0.40 0.20 0.10 0.00 0.30 0.00;]; 
Info1=zeros(num_Period,num_Object);

%各周期内传感器的融合?
for i=1:num_Period
Info1(i,:)=Info(1,:,i);
for j=1:num_Sensor-1
Info1(i,:)=DS_fusion(Info1(i,:),Info(j+1,:,i));
end
end
%周期之间的融合?
Result=Info1(1,:);
for i=1:num_Period-1
Result=DS_fusion(Result,Info1(i+1,:));
end
ec1=0.1; %融合决策判据?
ec2=0.1;
DS_out(Result,ec1,ec2);


function x=DS_fusion(x,y)
% 功能：融合x,y两行向量(经典Dempster-Shafer组合公式)
% x,y的格式形如[m1 m2 m3, ... , mk, m(全集), m(空集)]
% 要求m1 m2 m3 ...之间互相无交集
% m(全集)可不为0，表示不确定度
% m(空集)肯定是0
[nx,mx]=size(x);
if 1~=nx
    disp('x应为行向量');
    return;
end
[ny,my]=size(y);
if 1~=ny
    disp('y应为行向量');
    return;
end
if mx~=my
    disp('x,y列数应相等');
    return;
end
temp=0;
for i=1:mx-1
    
    if i==mx-1
        x(1,i)=x(1,i)*y(1,i);  %对全集的特殊处理
    else
        x(1,i)=x(1,i)*y(1,i)+x(1,i)*y(1,mx-1)+y(1,i)*x(1,mx-1);
    end
    temp=temp+x(1,i);
end
for i=1:mx-1
    x(1,i)=x(1,i)/temp;
end
x(1,mx)=0;
end