% ��UKF�Ľ��������˲��㷨--UPF
% ��UKF��������ֲ�
% �������˵����
%    Xiset����t-1ʱ�̵����Ӽ��ϣ�Z��tʱ�̵Ĺ۲�
%    Pin��ӦXiset���Ӽ��ϵķ���\
% �������˵����
%    Xo��upf�㷨���յĹ��ƽ��
%    Xoset��kʱ�̵����Ӽ��ϣ����ֵ����Xo
%    Pout��Xoset��Ӧ�ķ���
function [Xo,Xoset,Pout]=upf(Xiset,Z,t,Pin,N,R,Qukf,Rukf,g1,g2)
 
resamplingScheme=1;

 
Xukf=ones(1,N);     
Xset_pre=ones(1,N);  
Zpre=ones(1,N); 
for i=1:N
 
    [Xukf(i),Pout(i)]=ukf(Xiset(i),Z,Pin(i),Qukf,Rukf,t);
  
    Xset_pre(i) = Xukf(i) + sqrtm(Pout(i))*randn;
end

 
for i=1:N
 
    Zpre(i) = feval('hfun',Xset_pre(i),t);
 
    lik = inv(sqrt(R)) * exp(-0.5*inv(R)*((Z-Zpre(i))^(2)))+1e-99;
    prior = ((Xset_pre(i)-Xiset(i))^(g1-1)) * exp(-g2*(Xset_pre(i)-Xiset(i)));
    proposal = inv(sqrt(Pout(i))) * ...
        exp(-0.5*inv(Pout(i)) *((Xset_pre(i)-Xukf(i))^(2)));
    w(i) = lik*prior/proposal;
end
 
w = w./sum(w);

 
if resamplingScheme == 1
    outIndex = residualR(1:N,w');
elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w');
else
    outIndex = multinomialR(1:N,w');
end

 
Xoset = Xset_pre(outIndex);
Pout = Pout(outIndex);
 
Xo = mean(Xoset);