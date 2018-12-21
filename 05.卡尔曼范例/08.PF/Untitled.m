% ��������
N = 100;   %��������
Q = 5;      %��������
R = 5;      %��������
T = 20;     %����ʱ��
theta = pi/T;       %��ת�Ƕ�
distance = 80/T;    %ÿ���ߵľ���
WorldSize = 100;    %�����С
X = zeros(2, T);    %�洢ϵͳ״̬
Z = zeros(2, T);    %�洢ϵͳ�Ĺ۲�״̬
P = zeros(2, N);    %��������Ⱥ
PCenter = zeros(2, T);  %�������ӵ�����λ��
w = zeros(N, 1);         %ÿ�����ӵ�Ȩ��
err = zeros(1,T);     %���
X(:, 1) = [50; 20];     %��ʼϵͳ״̬
%wgn��m,n,p������һ��m��n�еĸ�˹�������ľ���p��dBWΪ��λָ�����������ǿ�ȡ�
Z(:, 1) = [50; 20] + wgn(2, 1, 10*log10(R));    %��ʼϵͳ�Ĺ۲�״̬

%��ʼ������Ⱥ
for i = 1 : N
    P(:, i) = [WorldSize*rand; WorldSize*rand];%��worldSize�������������N������
    dist = norm(P(:, i)-Z(:, 1));     %�����λ�����ľ��룬���ڹ�������ӵ�Ȩ��
    %���������Ѿ����������N�����ӣ����ڽ�������ʵ�Ĳ���ֵz���бȽϣ�Խ�ӽ���Ȩ��Խ�󣬻���˵��ֵԽСȨ��Խ��
    %�����Ȩ�ؼ����ǹ���p(z/x)�ķֲ������۲ⷽ�̵ķֲ�������۲����������˹�ֲ�����ôw��i��=p(z/x)
    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %��Ȩ��
end
PCenter(:, 1) = sum(P, 2) / N;      %�������ӵļ�������λ��

%%
err(1) = norm(X(:, 1) - PCenter(:, 1));     %���Ӽ���������ϵͳ��ʵ״̬�����
figure(1);
%����set(gca,'propertyname','propertyvalue'......)������Ե���ͼ�ε��������ԡ�
%propertyname��property value�ֱ�Ϊ�����������ƣ���������ĸ��Сд������Ӧ������ֵ��
set(gca,'FontSize',10);
hold on%ʹ��ǰ��ͼ�α��ֶ�����ˢ�£�׼�����ܴ˺󽫻��Ƶ�������
plot(X(1, 1), X(2, 1), 'r.', 'markersize',30)   %ϵͳ״̬λ��
%axis[xmin,xmax,ymin,ymax]�趨���귶Χ����������xmin<xmax������ͬ��
axis([0 100 0 100]);
plot(P(1, :), P(2, :), 'kx', 'markersize',5);   %�����������λ��
plot(PCenter(1, 1), PCenter(2, 1), 'b.', 'markersize',25); %�������ӵļ�������λ��
%ͼ�α�ע����legend(string1,string2,string3...),�ڵ�ǰͼ�����ͼ��

legend('True State', 'Particles', 'The Center of Particles');
%���ͼ�α�������title('string'),�ڵ�ǰ����ϵ�Ķ�����һ���ı���string,��Ϊ��ͼ�εı���
title('Initial State');
%ʹ��ǰ�ἰͼ�β��پ߱�����ˢ�µ�����
hold off

%%
%��ʼ�˶�
for k = 2 : T

    %ģ��һ�������˶���״̬
    %wgn��m,n,p������һ��m��n�еĸ�˹�������ľ���p��dBWΪ��λָ�����������ǿ�ȡ�
    X(:, k) = X(:, k-1) + distance * [(-cos(k * theta)); sin(k * theta)] + wgn(2, 1, 10*log10(Q));     %״̬����
    Z(:, k) = X(:, k) + wgn(2, 1, 10*log10(R));     %�۲ⷽ�� 

    %�����˲�
    %Ԥ��
    for i = 1 : N
        %��֮ǰ���ɵ����Ӵ��뵽״̬���̣�������һ��״̬��Ԥ��
        P(:, i) = P(:, i) + distance * [-cos(k * theta); sin(k * theta)] + wgn(2, 1, 10*log10(Q));
        dist = norm(P(:, i)-Z(:, k));     %�����λ�����ľ��룬���ڹ�������ӵ�Ȩ��
        %���������Ѿ����������N�����ӣ����ڽ�������ʵ�Ĳ���ֵz���бȽϣ�Խ�ӽ���Ȩ��Խ�󣬻���˵��ֵԽСȨ��Խ��
        %�����Ȩ�ؼ����ǹ���p(z/x)�ķֲ������۲ⷽ�̵ķֲ�������۲����������˹�ֲ�����ôw��i��=p(z/x)
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %��Ȩ��
    end
%��һ��Ȩ��
    wsum = sum(w);
    for i = 1 : N
        w(i) = w(i) / wsum;
    end

    %�ز��������£�
    for i = 1 : N
        wmax = 2 * max(w) * rand;  %��һ���ز�������,��ȡȨ�ص��е����ֵ����ȷֲ���0-1֮�����ĳ˻��ٳ���2��Ϊ�����ж�ѡ����Ȩ�����ӵ�����
        %randi()�������ɾ��ȷֲ���α�����������ΧΪimin--imax�����ûָ��imin����Ĭ��Ϊ1
        %r = randi([imin,imax],...)����һ����[imin,imax]��Χ�ڵ�α�������
        index = randi(N, 1);
        %�������������N�����ӵ��е�ĳ�����ӣ�Ȼ��ѡ������ӵ�Ȩֵ���ж��Ƿ��Ǵ������Ȩ�أ��ǵĻ��ͰѸ�����ѡ����
        while(wmax > w(index))
            %ʹvmax�ݼ�����Ȼ�Ļ��������ӵ�Ȩ���ǲ����ܴ���vmax��ֵ��
            wmax = wmax - w(index);
            index = index + 1;
            %�����ĳ����������ӿ�ʼ�ң�û�ҵ��Ӹ����ӿ�ʼ����������Ȩֵ�ܺʹ���vmax����ô�����´ӵ�һ�����ӿ�ʼ����
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     %�������index�еõ�������
    end
    %sum(X,2)��ʾ��X�������;�����sum(X),�Ǿ��ǰ������
    PCenter(:, k) = sum(P, 2) / N;      %�������ӵļ�������λ��

    %�������
    err(k) = norm(X(:, k) - PCenter(:, k));     %���Ӽ���������ϵͳ��ʵ״̬�����

    figure(2);
    set(gca,'FontSize',12);
    %clf; �������ͼ�ε����һ���ڻ�ͼ֮ǰ��
    clf;
    hold on
    plot(X(1, k), X(2, k), 'r.', 'markersize',50);  %ϵͳ״̬λ��
    axis([0 100 0 100]);
    plot(P(1, :), P(2, :), 'k.', 'markersize',5);   %��������λ��
    plot(PCenter(1, k), PCenter(2, k), 'b.', 'markersize',25); %�������ӵ�����λ��
    legend('True State', 'Particle', 'The Center of Particles');
    hold off
    pause(0.1);
end

%%
figure(3);
set(gca,'FontSize',12);
plot(X(1,:), X(2,:), 'r', Z(1,:), Z(2,:), 'g', PCenter(1,:), PCenter(2,:), 'b-');
axis([0 100 0 100]);
legend('True State', 'Measurement', 'Particle Filter');
xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);

%%
figure(4);
set(gca,'FontSize',12);
plot(err,'.-');
xlabel('t', 'FontSize', 20);
title('The err');