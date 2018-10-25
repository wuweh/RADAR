function [xe,xee,pk]=mykalmanadfun(A,U,C,Q,R,xe,y,p,Ea)

   xee=A*xe+U*Ea;
   p=A*p*A'+Q;
   l=p*C'*inv(C*p*C'+R);
   xe=xee+l*(y-C*xe);
   pk=inv(inv(p)+C'*inv(R)*C);