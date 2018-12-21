%%����������̵������˲�
clear all;close all;clc
% ��������
N = 200;   %��������
t = 20;                     % ����ʱ��
Ts = 0.1;                   % �������� 
len = fix(t/Ts);            % ���沽��
kx = .01;   ky = .05;       % ����ϵ��
g = 9.8;                    % ����
X(1,:) = [0, 50, 500, 0]; % ״̬ģ��ĳ�ֵ
dax = 2; day = 2;       % ϵͳ����

%������ʵ�켣
for k=2:len
    x = X(k-1,1); vx = X(k-1,2); y = X(k-1,3); vy = X(k-1,4); 
    x = x + vx*Ts;
    vx = vx + (-kx*vx^2+dax*randn)*Ts;
    y = y + vy*Ts;
    vy = vy + (ky*vy^2-g+day*randn)*Ts;
    X(k,:) = [x, vx, y, vy];
end

%��������
mrad = 0.001;
dr = 4; 
dafa = 10*mrad; % ��������
for k=1:len
    r = sqrt(X(k,1)^2+X(k,3)^2) + dr*randn(1,1);
    a = atan(X(k,1)/X(k,3)) + dafa*randn(1,1);
    Z(k,:) = [r, a];
end

V = [1,0.01,1,0.01];
% ��һ����˹�ֲ�����Ĳ�����ʼ������  
for i = 1:N    
    x_P(i,:) = X(1,:) + sqrt(V) * randn;    
end

PCenter(1,:) = mean(x_P); 

%��ʼ�˶�
for k = 2: len
    %�����˲�
    %Ԥ��
    %��ÿһ�����ӽ���һ��״̬����
    for i = 1 : N
        x = x_P(i,1); vx = x_P(i,2); y = x_P(i,3); vy = x_P(i,4);
        x = x + vx*Ts;
        vx = vx + (-kx*vx^2+dax*randn(1,1))*Ts;
        y = y + vy*Ts;
        vy = vy + (ky*vy^2-g+day*randn(1))*Ts;
        x_P_update(i,:) = [x, vx, y, vy];
        dist = sqrt((x_P_update(i,1)-Z(k,1)*sin(Z(k,2)))^2+(x_P_update(i,3)-Z(k,1)*cos(Z(k,2)))^2);       
        w(i) = (1/sqrt(dr)/sqrt(2*pi)) * exp(-(dist)^2/2/dr);   %��Ȩ��
    end
    
    %��һ��Ȩ��
    w = w./sum(w);

    %�ز���������һ��
%     for i = 1 : N
%         wmax = 2 * max(w) * rand;  %��һ���ز�������,��ȡȨ�ص��е����ֵ����ȷֲ���0-1֮�����ĳ˻��ٳ���2��Ϊ�����ж�ѡ����Ȩ�����ӵ�����
%         %randi()�������ɾ��ȷֲ���α�����������ΧΪimin--imax�����ûָ��imin����Ĭ��Ϊ1
%         %r = randi([imin,imax],...)����һ����[imin,imax]��Χ�ڵ�α�������
%         index = randi(N, 1);
%         %�������������N�����ӵ��е�ĳ�����ӣ�Ȼ��ѡ������ӵ�Ȩֵ���ж��Ƿ��Ǵ������Ȩ�أ��ǵĻ��ͰѸ�����ѡ����
%         while(wmax > w(index))
%             %ʹvmax�ݼ�����Ȼ�Ļ��������ӵ�Ȩ���ǲ����ܴ���vmax��ֵ��
%             wmax = wmax - w(index);
%             index = index + 1;
%             %�����ĳ����������ӿ�ʼ�ң�û�ҵ��Ӹ����ӿ�ʼ����������Ȩֵ�ܺʹ���vmax����ô�����´ӵ�һ�����ӿ�ʼ����
%             if index > N
%                 index = 1;
%             end          
%         end
%         x_P(i,:) = x_P_update(index,:);     %�������index�еõ�������
%     end
    
    %�ز�������������
    for i = 1 : N    
        x_P(i,:) = x_P_update(find(rand <= cumsum(w),1),:);   % ����Ȩ�ش�Ľ���õ����    
    end  

    PCenter(k,:) = mean(x_P); 
    %�������
    err(k) = sqrt((PCenter(k,1)-X(k,1))^2 + (PCenter(k,3)-X(k,3))^2);     %���Ӽ���������ϵͳ��ʵ״̬�����

%     figure(1);
%     set(gca,'FontSize',12);
    plot(X(k,1), X(k,3), 'r*', 'markersize',3);    hold on  %ϵͳ״̬λ��
    plot(x_P(:,1), x_P(:,3), 'k.', 'markersize',5);     hold on  %��������λ��
    plot(PCenter(k,1), PCenter(k,3), 'b*', 'markersize',3);    hold on ;grid on%�������ӵ�����λ��
    axis([0 300 0 500])
%     legend('True State', 'Particle', 'The Center of Particles');
     pause(0.001);
end

%��ͼ����
figure;
set(gca,'FontSize',12);
plot(X(:,1), X(:,3), 'b.-');hold on;
plot(Z(:,1).*sin(Z(:,2)), Z(:,1).*cos(Z(:,2)),'.');hold on;
plot(PCenter(:,1), PCenter(:,3), 'r.-');hold on;
grid on;axis([0 350 0 550])
legend('True State', 'Measurement', 'Particle Filter');
xlabel('x', 'FontSize', 10); ylabel('y', 'FontSize', 10);

figure;
set(gca,'FontSize',10);
plot(err,'.-');
xlabel('t', 'FontSize', 10);
title('The err');