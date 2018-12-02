function [A,Q,qa]=Starmodel(T,xa,a,xamax)
if xa>0
    qa=(xamax-xa)^2*(4-pi)/pi;
else
     qa=(xamax+xa)^2*(4-pi)/pi;
end
 
% A=[1 T (a*T-1+exp(-a*T))/a/a;0 1 (1-exp(-a*T))/a;0 0 exp(-a*T)];
temp = zeros(3,3);
A=[1 T (a*T-1+exp(-a*T))/a/a;0 1 (1-exp(-a*T))/a;0 0 exp(-a*T)];
% a1 = cat(2,A ,temp);
% a2 = cat(2,temp ,A);
% A = cat(1 , a1,a2);
% A=[1 T (a*T-1+exp(-a*T))/a/a 0 0 0;0 1 (1-exp(-a*T))/a 0 0 ;0 0 exp(-a*T) 0 0 0; 0 0 0 1 T (a*T-1+exp(-a*T))/a/a; 0 0 0 0 1 (1-exp(-a*T))/a;0 0 0 0 0 exp(-a*T)];

q11=(1-exp(-2*a*T)+2*a*T+2*(a^3)*(T^3)/3-2*a^2*T^2-4*a*T*exp(-a*T))/(a^4);
q12=(exp(-2*a*T)+1-2*exp(-a*T)+2*a*T*exp(-a*T)-2*a*T+a^2*T^2)/(a^3);
q13=(1-exp(-2*a*T)-2*a*T*exp(-a*T))/(a^2);
q22=(4*exp(-a*T)-3-exp(-2*a*T)+2*a*T)/(a^2);
q23=(exp(-2*a*T)+1-2*exp(-a*T))/a;
q33=(1-exp(-2*a*T));

temp = zeros(3,3);
Q=qa*[q11 q12 q13;q12 q22 q23;q13 q23 q33];
% a1 = cat(2,Q ,temp);
% a2 = cat(2,temp ,Q);
% Q = cat(1 , a1,a2);