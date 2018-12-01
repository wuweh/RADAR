%<Ŀ�����ǰ��������Ӧ�� 3.6 P72>
%SIR�����˲���Ӧ�ã��㷨���̲μ�����http://blog.csdn.net/heyijia0327/article/details/40899819    
clear all    
close all    
clc    
%initialize the variables    
x = 0.1; % initial actual state    
x_N = 1; % ϵͳ����������Э����  (������һά�ģ�������Ƿ���)    
x_R = 1; % ������Э����    
T = 75;  % ������75��    
N = 200; % ��������Խ��Ч��Խ�ã�������ҲԽ��    
    
%initilize our initial, prior particle distribution as a gaussian around    
%the true initial value    
    
V = 2; %��ʼ�ֲ��ķ���    
x_P = []; % ���� 

% ��һ����˹�ֲ�����Ĳ�����ʼ������    
for i = 1:N    
    x_P(i) = x + sqrt(V) * randn;    
end    
    
z_out = x^2 / 20 + sqrt(x_R) * randn;  %ʵ�ʲ���ֵ    
x_out = x;  %the actual output vector for measurement values.    
x_est = x; % time by time output of the particle filters estimate    
x_est_out = x_est; % the vector of particle filter estimates.    
    
for t = 1:T    
    x = 0.5*x + 25*x/(1 + x^2) + 8*cos(1.2*(t-1)) + sqrt(x_N)*randn;    
    z = x^2/20 + sqrt(x_R)*randn;    
    for i = 1:N    
        %������p(x(k)|x(k-1))�в���    
        x_P_update(i) = 0.5*x_P(i) + 25*x_P(i)/(1 + x_P(i)^2) + 8*cos(1.2*(t-1)) + sqrt(x_N)*randn;    
        %����������ӵ�ֵ��Ϊ���������Ȼȥ����Ȩ�����̵�    
        z_update(i) = x_P_update(i)^2/20;    
        %��ÿ�����Ӽ�����Ȩ�أ�����������������Ǹ�˹�ֲ������� w = p(y|x)��Ӧ����ļ��㹫ʽ    
        P_w(i) = (1/sqrt(2*pi*x_R)) * exp(-(z - z_update(i))^2/(2*x_R));    
    end    
    % ��һ��.    
    P_w = P_w./sum(P_w);    
      
    %Resampling����û���ò�����֮ǰ˵��histc����������Ŀ�ĺ�Ч����һ����    
    for i = 1 : N    
        x_P(i) = x_P_update(find(rand <= cumsum(P_w),1));   % ����Ȩ�ش�Ľ���õ����    
    end                                                     % find( ,1) ���ص�һ������ǰ�������������±�    
        
    %״̬���ƣ��ز����Ժ�ÿ�����ӵ�Ȩ�ض������1/N    
    x_est = mean(x_P);    
        
    % Save data in arrays for later plotting    
    x_out = [x_out x];    
    z_out = [z_out z];    
    x_est_out = [x_est_out x_est];      
end    
    
t = 0:T;       
plot(t, x_out, '.-b', t, x_est_out, '-.r','linewidth',3);    
set(gca,'FontSize',12); set(gcf,'Color','White');    
xlabel('time step'); ylabel('flight position');    
legend('True flight position', 'Particle filter estimate');    
