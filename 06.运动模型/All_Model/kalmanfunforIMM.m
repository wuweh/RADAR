function [xe,pk,v,S]=kalmanfunforIMM(A,C,Q,R,xe,y,p)
%This function is to calculate the estimation state and the real state.
   xe=A*xe;
   v=y-C*xe;
   p=A*p*A'+Q;
   S=C*p*C'+R;
   l=p*C'*inv(S);
   xe=xe+l*v;
   pk=(eye(size(p))-l*C)*p;
end
   
  