%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �����˲�һάϵͳ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Particle_For_UnlineOneDiv
clear all;close all;clc;
randn('seed',1); %Ϊ�˱�֤ÿ�����н��һ�£���������������ӵ�
%��ʼ����ز���
T=50;%��������
dt=1;%��������
Q=10;%������������
R=1;%������������
v=sqrt(R)*randn(T,1);%��������
w=sqrt(Q)*randn(T,1);%��������
numSamples=100;%������
ResampleStrategy=2;%=1Ϊ���������=2Ϊϵͳ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=0.1;%��ʼ״̬
%������ʵ״̬�͹۲�ֵ
X=zeros(T,1);%��ʵ״̬
Z=zeros(T,1);%����
X(1,1)=x0;%��ʵ״̬��ʼ��
Z(1,1)=(X(1,1)^2)./20+v(1,1);%�۲�ֵ��ʼ��
for k=2:T
	%״̬����
	X(k,1)=0.5*X(k-1,1)+2.5*X(k-1,1)/(1+X(k-1,1)^2)+8*cos(1.2*k)+w(k-1,1);
	%�۲ⷽ��
	Z(k,1)=(X(k,1).^2)./20+v(k,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����˲�����ʼ������Ҫ�������ڴ���˲�����״̬�����Ӽ��ϣ�Ȩ�ص�����
Xpf=zeros(numSamples,T);%�����˲�����״̬
Xparticles=zeros(numSamples,T);%���Ӽ���
Zpre_pf=zeros(numSamples,T);%�����˲��۲�Ԥ��ֵ
weight=zeros(numSamples,T);%Ȩ�س�ʼ��
%����״̬�͹۲�Ԥ��ĳ�ʼ������
Xpf(:,1)=x0+sqrt(Q)*randn(numSamples,1);
Zpre_pf(:,1)=Xpf(:,1).^2/20;
%������Ԥ�����
for k=2:T
	%��һ�������Ӽ��ϲ�������
	for i=1:numSamples
		QQ=Q;%���������˲���ͬ�������Q��Ҫ���������������һ��
		net=sqrt(QQ)*randn;%�����QQ���Կ��������İ뾶����ֵ�ɵ�
		Xparticles(i,k)=0.5.*Xpf(i,k-1)+2.5.*Xpf(i,k-1)./(1+Xpf(i,k-1).^2)+8*cos(1.2*k)+net;
	end
	%�ڶ����������Ӽ����е�ÿ�����ӣ���������Ҫ��Ȩֵ
	for i=1:numSamples
		Zpre_pf(i,k)=Xparticles(i,k)^2/20;
		weight(i,k)=exp(-.5*R^(-1)*(Z(k,1)-Zpre_pf(i,k))^2);%ʡ���˳�����
	end
	weight(:,k)=weight(:,k)./sum(weight(:,k));%��һ��Ȩֵ
	%������������Ȩֵ��С�����Ӽ����ز�����Ȩֵ���Ϻ����Ӽ�����һһ��Ӧ��
	%ѡ���������
	if ResampleStrategy==1
		outIndex = randomR(weight(:,k));
	elseif ResampleStrategy==2
		outIndex = systematicR(weight(:,k)');
	elseif ResampleStrategy==3
		outIndex = multinomialR(weight(:,k));
	elseif ResampleStrategy==4
		outIndex = residualR(weight(:,k)');
	end
	%���Ĳ��������ز����õ���������ȥ��ѡ��Ӧ�����ӣ��ع��ļ��ϱ����˲����״̬����
	%�����״̬�������ֵ���������յ�Ŀ��״̬��
	Xpf(:,k)=Xparticles(outIndex,k);
end
%��������ֵ���ơ���������Ƽ����Ʒ���
Xmean_pf=mean(Xpf);%�����ֵ���ƣ�������ĵ��Ĳ���Ҳ�������˲����Ƶ�����״̬
bins=20;
Xmap_pf=zeros(T,1);
for k=1:T
	[p,pos]=hist(Xpf(:,k,1),bins);
	map=find(p==max(p));
	Xmap_pf(k,1)=pos(map(1));%���������
end
for k=1:T
	Xstd_pf(1,k)=std(Xpf(:,k)-X(k,1));%��������׼�����
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ͼ
figure();clf;%���������Ͳ�������ͼ
subplot(221);
plot(v);%��������
xlabel('ʱ��');ylabel('��������');
subplot(222);
plot(w);%��������
xlabel('ʱ��');ylabel('��������');
subplot(223);
plot(X);%��ʵ״̬
xlabel('ʱ��');ylabel('״̬X');
subplot(224);
plot(Z);%�۲�ֵ
xlabel('ʱ��');ylabel('�۲�Z');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
k=1:dt:T;
plot(k,X,'b',k,Xmean_pf,'r',k,Xmap_pf,'g');%ע��Xmean_pf���������˲����
legend('ϵͳ��ʵ״ֵ̬','�����ֵ����','��������ʹ���');
xlabel('ʱ��');ylabel('״̬����');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
subplot(121);
plot(Xmean_pf,X,'+');%�����˲�����ֵ����ʵ״ֵ̬���1:1��ϵ�����ԳƷֲ�
xlabel('�����ֵ����');ylabel('��ֵ');
hold on;
c=-25:1:25;
plot(c,c,'r');%����ɫ�ĶԳ���y=x
hold off;
subplot(122);%���������ֵ����ʵ״ֵ̬���1:1��ϵ�����ԳƷֲ�
plot(Xmap_pf,X,'+');
xlabel('Map����');ylabel('��ֵ');
hold on;
c=-25:25;
plot(c,c,'r');%����ɫ�ĶԳ���y=x
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ֱ��ͼ����ͼ����Ϊ�˿����Ӽ��ĺ����ܶ�
domain=zeros(numSamples,1);
range=zeros(numSamples,1);
bins=10;
support=[-20:1:20];
figure();
hold on;%ֱ��ͼ
xlabel('ʱ��');ylabel('�����ռ�');
vect=[0 1];
caxis(vect);
for k=1:T
	%ֱ��ͼ��ӳ�˲�������Ӽ��ϵķֲ����
	[range,domain]=hist(Xpf(:,k),support);
	%����waterfall��������ֱ��ͼ�ֲ������ݻ�����
	waterfall(domain,k,range);
end
axis([-20 20 0 T 0 100]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
xlabel('�����ռ�');ylabel('�����ܶ�');
k=30;%k=?��ʾҪ�鿴�ڼ���ʱ�̵����ӷֲ�����ʵ״ֵ̬���ص���ϵ
[range,domain]=hist(Xpf(:,k),support);
plot(domain,range);
%��ʵ״̬�������ռ��е�λ�ã���һ����ɫֱ�߱�ʾ
XXX=[X(k,1),X(k,1)];
YYY=[0,max(range)+10];
line(XXX,YYY,'Color','r');
axis([min(domain) max(domain) 0 max(range)+10]);
figure();
k=1:dt:T;
plot(k,Xstd_pf,'-');
xlabel('ʱ��');ylabel('״̬��������׼��');
axis([0,T,0,10]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������ܣ�ʵ������ز����㷨
%���������weightΪԭʼ���ݶ�Ӧ��Ȩ�ش�С
%���������outIndex�Ǹ���weight��inIndexɸѡ�͸��ƽ��
function outIndex=randomR(weight)
%������ݵĳ���
L=length(weight);
%��ʼ������������������������������������
outIndex=zeros(1,L);
%��һ��������[0,1]�Ͼ��ȷֲ���������飬����������
u=unifrnd(0,1,1,L);
u=sorf(u);
%u=(1:L)/L%�������ȫ����
%�ڶ�������������Ȩ�ػ��ۺ���cdf
cdf=cumsum(weight);
%�����������ļ���
i=1;
for j=1:L
	%�˴��Ļ���ԭ���ǣ�u�Ǿ��ȵģ���Ȼ��Ȩֵ��ĵط�
	%�и�����������������䣬��˻ᱻ��θ���
	while(i<=L)&(u(i)<=cdf(j))
		%����Ȩֵ�������
		outIndex(i)=j;
		%����������һ������������������ĸ�����
		i=i+1;
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����ʽ�ز����Ӻ���
% ���������weightΪԭʼ���ݶ�Ӧ��Ȩ�ش�С
% ���������outIndex�Ǹ���weightɸѡ�͸��ƽ��
function outIndex = multinomialR(weight)
%��ȡ���ݳ���
Col=length(weight);
N_babies= zeros(1,Col);

%��������Ȩ���ۼƺ���cdf 
cdf= cumsum(weight);
 %����[0,1]���ȷֲ��������
u=rand(1,Col);

%��u^(j^-1)�η� 
uu=u.^(1./(Col:-1:1));
 %���A��һ��������cumprod(A)������һ������A��Ԫ�ػ������˵Ľ��������
 %Ԫ�ظ�����ԭ������ͬ
ArrayTemp=cumprod(uu);
 %fliplr(X)ʹ����X�ش�ֱ�����ҷ�ת
u = fliplr(ArrayTemp);
j=1;
for i=1:Col
    %�˴��������������
    while (u(i)>cdf(j))
        j=j+1;
    end
    N_babies(j)=N_babies(j)+1;
end
index=1;
for i=1:Col
    if (N_babies(i)>0)
        for j=index:index+N_babies(i)-1
            outIndex(j) = i;
        end
    end
    index= index+N_babies(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��������˵�����в��ز�������
% ���������һ��Ȩ��weight����
% ���������Ϊ��Ȩ���ز����������outIndex
function outIndex = residualR(weight)
N= length(weight);
N_babies= zeros(1,N);
q_res = N.*weight;
N_babies = fix(q_res);
N_res=N-sum(N_babies);
if (N_res~=0)
    q_res=(q_res-N_babies)/N_res;
    cumDist= cumsum(q_res);
    u = fliplr(cumprod(rand(1,N_res).^(1./(N_res:-1:1))));
    j=1;
    for i=1:N_res
        while (u(1,i)>cumDist(1,j))
            j=j+1;
        end
        N_babies(1,j)=N_babies(1,j)+1;
    end
end
index=1;
for i=1:N
    if (N_babies(1,i)>0)
        for j=index:index+N_babies(1,i)-1
            outIndex(j) = i;
        end
    end
    index= index+N_babies(1,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ϵͳ�ز����Ӻ���
% ���������weightΪԭʼ���ݶ�Ӧ��Ȩ�ش�С
% ���������outIndex�Ǹ���weightɸѡ�͸��ƽ��
function outIndex = systematicR(weight)
N=length(weight);
N_children=zeros(1,N);
label=zeros(1,N);
label=1:1:N;
s=1/N;
auxw=0;
auxl=0;
li=0;
T=s*rand(1);
j=1;
Q=0;
i=0;
u=rand(1,N);
while (T<1)
    if (Q>T)
        T=T+s;
        N_children(1,li)=N_children(1,li)+1;
    else
        i=fix((N-j+1)*u(1,j))+j;
        auxw=weight(1,i);
        li=label(1,i);
        Q=Q+auxw;
        weight(1,i)=weight(1,j);
        label(1,i)=label(1,j);
        j=j+1;
    end
end
index=1;
for i=1:N
    if (N_children(1,i)>0)
        for j=index:index+N_children(1,i)-1
            outIndex(j) = i;
        end
    end
    index= index+N_children(1,i);
end