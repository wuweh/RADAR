% author :  Perry.Li  @USTC
% function: simulating the process of EKF
% date:     04/28/2015
% 
N = 50;         %��������N��ʱ�� 
n=3;            %״̬ά��
q=0.1;          %���̱�׼��
r=0.2;          %������׼��
Q=q^2*eye(n);   %���̷���
R=r^2;          %����ֵ�ķ��� 
f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  %״̬����
h=@(x)[x(1);x(2);x(3)];                   %��������
s=[0;0;1];                                %��ʼ״̬
%��ʼ��״̬
x=s+q*randn(3,1);                         
P = eye(n);                               
xV = zeros(n,N);          
sV = zeros(n,N);         
zV = zeros(n,N);
for k=1:N
  z = h(s) + r*randn;                     
  sV(:,k)= s;                             %ʵ��״̬
  zV(:,k)  = z;                           %״̬����ֵ
  [x1,A]=jaccsd(f,x); %����f���ſɱȾ�������x1��Ӧ�ƽ�ʽline2
  P=A*P*A'+Q;         %���̷���Ԥ�⣬��Ӧline3
  [z1,H]=jaccsd(h,x1); %����h���ſɱȾ���
  K=P*H'*inv(H*P*H'+R); %���������棬��Ӧline4
  x=x1+K*(z-z1);        %״̬EKF����ֵ����Ӧline5
  P=P-K*H*P;            %EKF�����Ӧline6
  xV(:,k) = x;          %save
  s = f(s) + q*randn(3,1);  %update process 
end
for k=1:3
  FontSize=14;
  LineWidth=1;
  figure();
  plot(sV(k,:),'g-'); %������ʵֵ
  hold on;
  plot(xV(k,:),'b-','LineWidth',LineWidth) %�������Ź���ֵ
  hold on;
  plot(zV(k,:),'k+'); %����״̬����ֵ
  hold on;
  legend('��ʵ״̬', 'EKF���Ź��ƹ���ֵ','״̬����ֵ');
  xl=xlabel('ʱ��(����)');
  t=['״̬ ',num2str(k)] ;
  yl=ylabel(t);
  set(xl,'fontsize',FontSize);
  set(yl,'fontsize',FontSize);
  hold off;
  set(gca,'FontSize',FontSize);
end