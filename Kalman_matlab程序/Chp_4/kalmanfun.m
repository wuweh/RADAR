function [xe,pk,p1]=kalmanfun(A,C,Q,R,xe,z,p)
%This function is to calculate the estimation state by Kalman filter.
   xe=A*xe;
   p1=A*p*A'+Q;
   K=p1*C'*inv(C*p1*C'+R);
   xe=xe+K*(z-C*xe);
   pk=(eye(size(p1))-l*C)*p1;