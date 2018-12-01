% ��EKF�Ľ��������˲��㷨--EPF
% ��EKF��������ֲ�
% �������˵����
%    Xiset����t-1ʱ�̵����Ӽ��ϣ�Z��tʱ�̵Ĺ۲�
%    Pin��ӦXiset���Ӽ��ϵķ���\
% �������˵����
%    Xo��epf�㷨���յĹ��ƽ��
%    Xoset��kʱ�̵����Ӽ��ϣ����ֵ����Xo
%    Pout��Xoset��Ӧ�ķ���
function [Xo,Xoset,Pout]=epf(Xiset,Z,t,Pin,N,R,Qekf,Rekf,g1,g2)
  
resamplingScheme=1;
Zpre=ones(1,N);     
Xsetpre=ones(1,N);   
w = ones(1,N);      

Pout=ones(1,N);      
Xekf=ones(1,N);     
Xekf_pre=ones(1,N); 

 
for i=1:N
    [Xekf(i),Pout(i)]=ekf(Xiset(i),Z,Pin(i),t,Qekf,Rekf);
    Xsetpre(i)=Xekf(i)+sqrtm(Pout(i))*randn;
end

 
for i=1:N
    Zpre(i) = feval('hfun',Xsetpre(i),t);
    lik = inv(sqrt(R)) * exp(-0.5*inv(R)*((Z-Zpre(i))^(2)))+1e-99;
    prior = ((Xsetpre(i)-Xiset(i))^(g1-1)) * exp(-g2*(Xsetpre(i)-Xiset(i)));
    proposal = inv(sqrt(Pout(i))) * ...
        exp(-0.5*inv(Pout(i)) *((Xsetpre(i)-Xekf(i))^(2)));
    w(i) = lik*prior/proposal;
end
 
w= w./sum(w);

 
if resamplingScheme == 1
    outIndex = residualR(1:N,w');
elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w');
else
    outIndex = multinomialR(1:N,w');
end
 
Xoset = Xsetpre(outIndex);
Pout = Pout(outIndex);
 
Xo = mean(Xoset);