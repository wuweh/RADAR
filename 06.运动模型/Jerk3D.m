%��άJerkģ��------
%{
    
%}

clear all;
clc;
%�������,��ʼ��
T=0.02;         %����s
t0=T;           %��ʼʱ��
t_end=10;    %��ֹʱ��s
t=t0:T:t_end;      %����ʱ���

N=length(t);  %����ʱ���

x=zeros(1,N);     %x�������ʼ��
y=zeros(1,N);     %y�������ʼ��
z=zeros(1,N);     %z�������ʼ��
x1=zeros(1,N);     %x���ٶȳ�ʼ��
y1=zeros(1,N);     %y���ٶȳ�ʼ��
z1=zeros(1,N);     %z���ٶȳ�ʼ��
x2=zeros(1,N);     %x����ٶȳ�ʼ��
y2=zeros(1,N);     %y����ٶȳ�ʼ��
z2=zeros(1,N);     %z����ٶȳ�ʼ��
x3=zeros(1,N);     %x��Ӽ��ٶȳ�ʼ��
y3=zeros(1,N);     %y��Ӽ��ٶȳ�ʼ��
z3=zeros(1,N);     %z��Ӽ��ٶȳ�ʼ��

xdelta_R=20;       %�״�x���������m
ydelta_R=30;       %�״�y���������m
zdelta_R=25;       %�״�z���������m

alpha=0.1;
jerk_x=-20;
g=-9.8;
X=[x;  x1; x2;  x3;    y;  y1;  y2;  y3;    z;  z1;  z2; z3];       %Ŀ��״̬����
%����ֵ
X(:,1)=[3000;  -200;  0;  jerk_x;    4000;  200;  0;  0;    10000;  0;  g; 0];
%{
x(1)=3000;          %x��ʼλ��m
y(1)=4000;          %y��ʼλ��m
vx(1)=200;          %x��ʼ�ٶ�m/s
vy(1)=200;          %y��ʼ�ٶ�m/s
vz(1)=0;            %z��ʼ�ٶ�m/s
%}

p1=(2 - 2*alpha*T + alpha^2*T^2 - 2*exp(-alpha*T))/(2*alpha^3);
q1=(exp(-alpha*T) - 1 + alpha*T)/alpha^2;
r1=(1 - exp(-alpha*T))/alpha;
s1=exp(-alpha*T);

%Ŀ���״̬����
F4=[ 1, T, T^2/2, p1;
     0, 1, T,     q1;
     0, 0, 1,     r1;
     0, 0, 0,     s1 ];
 
F=[F4,              zeros(4,4),     zeros(4,4);
   zeros(4,4),      F4,             zeros(4,4);
   zeros(4,4),      zeros(4,4),     F4];

 H=[ 1, 0, 0, 0,  0, 0, 0,  0,  0, 0, 0, 0;
        0, 0, 0,  0, 1, 0, 0, 0,   0, 0, 0, 0;
        0, 0, 0, 0,  0, 0, 0, 0,   1, 0, 0, 0  ];
Y=zeros(size(H,1),N);   %��ʼ������״̬.
r=0.05;
for k=2:N 
    X(:,k)=F*X(:,k-1);
    X(11,k)=g;
    Y(:,k)=H*X(:,k-1)+r*randn(size(H,1),1);
end
    
x=X(1,:);
y=X(5,:);
z=X(9,:);

for i=1:1:N
    R(i)=sqrt( X(1,i)^2 + X(5,i)^2 + X(9,i)^2 );
    v(i)=sqrt( X(2,i)^2 + X(6,i)^2 + X(10,i)^2 );
    a(i)=sqrt( X(3,i)^2 + X(7,i)^2 + X(11,i)^2 );
end

figure(1)
plot3(x,y,z);
xlabel('x��:m'); zlabel('y��:m');grid on;  

figure(2)
subplot(311),plot(t,R,'-b');
legend('R'),title('R����ͼ'); xlabel('ʱ�䣺s'); ylabel('R:m');
subplot(312),plot(t,v);
legend('v'),title('�ٶ�����'); xlabel('����ʱ���'); ylabel('�ٶ�:m/s');
subplot(313),plot(t,a);
legend('a'),title('���ٶ�����'); xlabel('t:s'); ylabel('���ٶ�:m/s^2');
