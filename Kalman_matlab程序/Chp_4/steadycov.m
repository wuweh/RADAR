function [pp1,pp,KK]=steadycov(A,C,Q,R,I,p0)
p=p0*ones(size(A));
pp1=[];pp=[];KK=[];
for i=1:I
p1=A*p*A'+Q;
K=p1*C'*inv(C'*p1*C+R);
p=(eye(size(A))-K*C)*p1;
pp1=[pp1 diag(p1)];
pp=[pp diag(p)];
KK=[KK K];
end