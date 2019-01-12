clear all
close all
clc

 %%%%%%%%%% ?D_S�㷨���� ? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%?
num_Sensor=3;  %��������Ŀ?
num_Object=6;  %ʶ�������Ŀ ?�ĸ�Ŀ��(�����״�: ?��� ��׼ ���� Ԥ��)+ȫ��+�ռ�?
num_Period=3;  %����������?

Info=zeros(num_Sensor,num_Object,num_Period);  % һ���о�����Ҫ����Ϣ?
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

%�������ڴ��������ں�?
for i=1:num_Period
Info1(i,:)=Info(1,:,i);
for j=1:num_Sensor-1
Info1(i,:)=DS_fusion(Info1(i,:),Info(j+1,:,i));
end
end
%����֮����ں�?
Result=Info1(1,:);
for i=1:num_Period-1
Result=DS_fusion(Result,Info1(i+1,:));
end
ec1=0.1; %�ںϾ����о�?
ec2=0.1;
DS_out(Result,ec1,ec2);


function x=DS_fusion(x,y)
% ���ܣ��ں�x,y��������(����Dempster-Shafer��Ϲ�ʽ)
% x,y�ĸ�ʽ����[m1 m2 m3, ... , mk, m(ȫ��), m(�ռ�)]
% Ҫ��m1 m2 m3 ...֮�以���޽���
% m(ȫ��)�ɲ�Ϊ0����ʾ��ȷ����
% m(�ռ�)�϶���0
[nx,mx]=size(x);
if 1~=nx
    disp('xӦΪ������');
    return;
end
[ny,my]=size(y);
if 1~=ny
    disp('yӦΪ������');
    return;
end
if mx~=my
    disp('x,y����Ӧ���');
    return;
end
temp=0;
for i=1:mx-1
    
    if i==mx-1
        x(1,i)=x(1,i)*y(1,i);  %��ȫ�������⴦��
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