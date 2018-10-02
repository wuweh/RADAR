%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ����˵���� UKF��Ŀ������е�Ӧ��
%  ����˵���� 1��״̬6ά��x�����λ�á��ٶȡ����ٶȣ�
%                y�����λ�á��ٶȡ����ٶȣ�
%             2���۲���ϢΪ����ͽǶȣ�
%  ��ϸԭ����ܼ�����ע����ο���
%  ���������˲�ԭ��Ӧ��-MATLAB���桷�����ӹ�ҵ�����磬��Сƽ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ukf_for_track_6_div_system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=6;
t=0.05;
%��������Э����
Q=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 0.01 0 0 0;
    0 0 0 0.01 0 0;
    0 0 0 0 0.0001 0;
    0 0 0 0 0 0.0001];
%��������Э����
R = [100 0;
    0 0.001^2];
%״̬����
%x1Ϊx��λ�ã�x2Ϊy��λ�ã�x3��x4�ֱ���x��y����ٶȣ�x5��x6�ֱ��Ǽ��ٶ�
f=@(x)[ x(1)+t*x(3)+0.5*t^2*x(5);
        x(2)+t*x(4)+0.5*t^2*x(6); 
        x(3)+t*x(5);
        x(4)+t*x(6);
        x(5);
        x(6)];
%�۲ⷽ��
h=@(x)[ sqrt(x(1)^2+x(2)^2);
        atan(x(2)/x(1))];

s1=[92000;82000;-130;-230;0;0];%ʵ��λ�ó�ʼ��
s2=[91000;82000;-130;-230;0;0];%ʵ��λ�ó�ʼ��
s3=[90000;82000;-130;-230;0;0];%ʵ��λ�ó�ʼ��
s4=[89000;82000;-130;-230;0;0];%ʵ��λ�ó�ʼ��
s5=[90550;90930;-250;190;32;-23];%ʵ��λ�ó�ʼ��
s6=[90120;90600;-139;260;47;-24];%ʵ��λ�ó�ʼ��
s7=[90990;91390;-200;320;53;-25];%ʵ��λ�ó�ʼ��



x0=s1+sqrtm(Q)*randn(n,1);%����λ�ó�ʼ��״̬

%��ʼ��Э����
P0 =[100 0 0 0 0 0;
    0 100 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 0.0001 0;
    0 0 0 0 0 0.0001];

%�ܷ��沽��
N=50;

%ukf�˲�״̬��ʼ��
Xukf = zeros(n,N);
X1 = zeros(n,N);
X2 = zeros(n,N);
X3 = zeros(n,N);
X4 = zeros(n,N);
X5 = zeros(n,N);
X6 = zeros(n,N);
X7 = zeros(n,N);
Z = zeros(2,N);

figure
hold on;
for i=1:N
    X1(:,i)= f(s1)+sqrtm(Q)*randn(6,1); %ģ�⣬����Ŀ���˶���ʵ�켣
    s1 = X1(:,i);
     s1_1(i) = sqrt(X1(1,i)^2+X1(2,i)^2);
    
    X2(:,i)= f(s2)+sqrtm(Q)*randn(6,1); %ģ�⣬����Ŀ���˶���ʵ�켣
    s2 = X2(:,i);
     s2_1(i) = sqrt(X2(1,i)^2+X2(2,i)^2);
    
    X3(:,i)= f(s3)+sqrtm(Q)*randn(6,1); %ģ�⣬����Ŀ���˶���ʵ�켣
    s3 = X3(:,i);
     s3_1(i) = sqrt(X3(1,i)^2+X3(2,i)^2);
    
    X4(:,i)= f(s4)+sqrtm(Q)*randn(6,1); %ģ�⣬����Ŀ���˶���ʵ�켣
    s4 = X4(:,i);
     s4_1(i) = sqrt(X4(1,i)^2+X4(2,i)^2);
    
    X5(:,i)= f(s5)+sqrtm(Q)*randn(6,1); %ģ�⣬����Ŀ���˶���ʵ�켣
    s5 = X5(:,i);
     s5_1(i) = sqrt(X5(1,i)^2+X5(2,i)^2);
    
    X6(:,i)= f(s6)+sqrtm(Q)*randn(6,1); %ģ�⣬����Ŀ���˶���ʵ�켣
    s6 = X6(:,i);
     s6_1(i) = sqrt(X6(1,i)^2+X6(2,i)^2);
    
    X7(:,i)= f(s7)+sqrtm(Q)*randn(6,1); %ģ�⣬����Ŀ���˶���ʵ�켣
    s7 = X7(:,i);
     s7_1(i) = sqrt(X7(1,i)^2+X7(2,i)^2); 
     
    plot(i,s1_1(i),'.g');
    plot(i,s2_1(i),'.r');
    plot(i,s3_1(i),'.y');
    plot(i,s4_1(i),'.b');
    plot(i,s5_1(i),'*g');
    plot(i,s6_1(i),'*r');
    plot(i,s7_1(i),'*y');
end
figure;
hold on;
plot(X1(1,:),X1(2,:),'.g');
plot(X2(1,:),X2(2,:),'.r');
plot(X3(1,:),X3(2,:),'.y');
plot(X4(1,:),X4(2,:),'.b');
plot(X5(1,:),X5(2,:),'*g');
plot(X6(1,:),X6(2,:),'*r');
plot(X7(1,:),X7(2,:),'*y');

ux=x0;
for k=1:N
    %h(X(:,k))ָ����ͨ�������õ�Ŀ���λ����Ϣ������ͽǶ�
    Z(:,k)= h(X1(:,k)) + sqrtm(R)*randn(2,1);        %����ֵ��������
    [Xukf(:,k), P0] = ukf(f,ux,P0,h,Z(:,k),Q,R);    %����ukf�㷨
    ux=Xukf(:,k);
end

%����������
for k=1:N
    RMS(k)=sqrt((X1(1,k)-Xukf(1,k))^2+(X1(2,k)-Xukf(2,k))^2 );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
t=1:N;
hold on;box on;
plot(X1(1,t),X1(2,t), 'k-')  %��ʵĿ��xy��Ϣ
plot(Z(1,t).*cos(Z(2,t)),Z(1,t).*sin(Z(2,t)),'-b.')%����Ŀ����Ϣ
plot(Xukf(1,t),Xukf(2,t),'-r.')%ukfĿ����Ϣ
legend('ʵ��ֵ','����ֵ','ukf����ֵ');
xlabel('x����λ��/��');
ylabel('y����λ��/��');

figure
box on;
plot(RMS,'-ko','MarkerFace','r')
xlabel('t/��')
ylabel('ƫ��/��')

%%
%ffun��״̬ת�ƾ���
%X������״̬
%P��Э����
% hfun:�۲ⷽ��
% Z:��ǰ�۲�״̬
% Q:��������Э����
% R:��������Э����
function [X,P]=ukf(ffun,X,P,hfun,Z,Q,R)
L=numel(X);  %���X��Ԫ�ظ���
m=numel(Z);  %���Z��Ԫ�ظ���
%������������
alpha=1e-2;
ki=0;
beta=2;

lambda=alpha^2*(L+ki)-L; %��ʽ1
c=L+lambda;              %�м����
Wm=[lambda/c 0.5/c+zeros(1,2*L)]; %��ʽ2
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);%��ʽ3

c=sqrt(c);
Xsigmaset=sigmas(X,P,c); 
[X1means,X1,P1,X2]=ut(ffun,Xsigmaset,Wm,Wc,L,Q);   
[Zpre,Z1,Pzz,Z2]=ut(hfun,X1,Wm,Wc,m,R);

Pxz=X2*diag(Wc)*Z2';
K=Pxz/(Pzz); %�˲�����
X=X1means+K*(Z-Zpre);   %����״̬��ֵ
P=P1-K*Pxz';            %����״̬Э����

%%
function [Xmeans,Xsigma_pre,P,Xdiv]=ut(fun,Xsigma,Wm,Wc,n,COV)
LL=size(Xsigma,2);
Xmeans=zeros(n,1);
Xsigma_pre=zeros(n,LL);
for k=1:LL
    Xsigma_pre(:,k)=fun(Xsigma(:,k));  %״̬��sigma������Ա任
    Xmeans=Xmeans+Wm(k)*Xsigma_pre(:,k); %״̬��ֵ
end
Xdiv=Xsigma_pre-Xmeans(:,ones(1,LL));
P=Xdiv*diag(Wc)*Xdiv'+COV;              %״̬Э����

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sigma points around reference point
% ����2n+1��sigma��
% Inputs:
% x: reference point
% P: covariance
% c: coefficient
% Output:
% X: Sigma points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xset=sigmas(X,P,c)
A = c*chol(P)';
Y = X(:,ones(1,numel(X)));
Xset = [X Y+A Y-A];
