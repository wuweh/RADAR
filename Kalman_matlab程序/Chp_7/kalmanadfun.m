function [xe,pk]=kalmanadfun(A1,A,U,C,Q,R,xe,y,p)
%This function is to calculate the estimation state and the real state.


   xe=A1*xe;
   p=A*p*A'+Q;
   l=p*C'*inv(C*p*C'+R);
   xe=xe+l*(y-C*xe);
   pk=(eye(size(p))-l*C)*p;
   
  