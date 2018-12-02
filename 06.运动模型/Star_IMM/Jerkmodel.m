function [A,Q]=Jerkmodel(T,alpha)
%状态噪声协方差矩阵------------------------------------------------------------
Q11=(alpha^5*T^5/10-alpha^4*T^4/2+4*alpha^3*T^3/3-2*alpha^2*T^2+2*alpha*T-3 ...
    +4*exp(-alpha*T)+2*alpha^2*T^2*exp(-alpha*T)-exp(-2*alpha*T))/(2*alpha^7);
Q12=(1-2*alpha*T+2*alpha^2*T^2-alpha^3*T^3+alpha^4*T^4/4+exp(-2*alpha*T)+ ...
    2*alpha*T*exp(-alpha*T)-2*exp(-alpha*T)-alpha^2*T^2*exp(-alpha*T))/(2*alpha^6);%与Q21相同,以下同理
Q13=(2*alpha*T-alpha^2*T^2+alpha^3*T^3/3-3-exp(-2*-alpha*T)+4*exp(-alpha*T) ...
    +alpha^2*T^2*exp(-alpha*T))/(2*alpha^5);
Q14=(1+exp(-2*alpha*T)-2*exp(-alpha*T)-alpha^2*T^2*exp(-alpha*T))/(2*alpha^4);
Q21=Q12;
Q22=(1-exp(-2*alpha*T)+2*alpha^3*T^3/3+2*alpha*T-2*alpha^2*T^2-4*alpha*T*exp(-alpha*T))/(2*alpha^5);
Q23=(1+alpha^2*T^2-2*alpha*T+2*alpha*T*exp(-alpha*T)+exp(-2*alpha*T)-2*exp(-alpha*T))/(2*alpha^4);
Q24=(1-2*alpha*T*exp(-2*alpha*T)-exp(-2*alpha*T))/(2*alpha^3);
Q31=Q13;
Q32=Q23;
Q33=(4*exp(-alpha*T)-exp(-2*alpha*T)+2*alpha*T-3)/(2*alpha^3);
Q34=(1-2*exp(-alpha*T)+exp(-2*alpha*T))/(2*alpha^2);
Q41=Q14;
Q42=Q24;
Q43=Q34;
Q44=(1-exp(-2*alpha*T))/(2*alpha);
 
Q=[ Q11,Q12,Q13,Q14;
    Q21,Q22,Q23,Q24;
    Q31,Q32,Q33,Q34;
    Q41,Q42,Q43,Q44];

%转移矩阵----------------------------------------------------------------------
PT=(2-2*alpha*T+alpha^2*T^2-2*exp(-alpha*T))/(2*alpha^3);
QT=(exp(-alpha*T)-1+alpha*T)/alpha^2;
RT=(1-exp(-alpha*T))/alpha;
ST=exp(-alpha*T);

A=[1,T,T^2/2,PT;
    0,1,T,QT;
    0,0,1,RT;
    0,0,0,ST];

end