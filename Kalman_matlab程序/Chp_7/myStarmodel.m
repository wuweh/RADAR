function [A,Q,U]=myStarmodel(T,a,qq)
 
A=[1 T (a*T-1+exp(-a*T))/a/a;0 1 (1-exp(-a*T))/a;0 0 exp(-a*T)];
 
q11=(1-exp(-2*a*T)+2*a*T+2*(a^3)*(T^3)/3-2*a^2*T^2-4*a*T*exp(-a*T))/2/(a^5);
q12=(exp(-2*a*T)+1-2*exp(-a*T)+2*a*T*exp(-a*T)-2*a*T+a^2*T^2)/2/(a^4);
q13=(1-exp(-2*a*T)-2*a*T*exp(-a*T))/2/(a^3);
q22=(4*exp(-a*T)-3-exp(-2*a*T)+2*a*T)/2/(a^3);
q23=(exp(-2*a*T)+1-2*exp(-a*T))/2/(a^2);
q33=(1-exp(-2*a*T))/2/a;
 
Q=qq*[q11 q12 q13;q12 q22 q23;q13 q23 q33];
 
U=[(-T+T^2/2+(1-exp(-a*T))/a)/a;T-(1-exp(-a*T))/a;1-exp(-a*T)];