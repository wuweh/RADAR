function [Xekf,Pout]=ekf(Xin,Z,Pin,t,Qekf,Rekf)

Xpre=feval('ffun',Xin,t);
Jx=0.5;
 
Pekfpre = Qekf + Jx*Pin*Jx';
Zekfpre= feval('hfun',Xpre,t);
 
if t<=30
    Jy = 2*0.2*Xpre;
else
    Jy = 0.5;
end
 
M = Rekf + Jy*Pekfpre*Jy';
K = Pekfpre*Jy'*inv(M);
Xekf=Xpre+K*(Z-Zekfpre);
Pout = Pekfpre - K*Jy*Pekfpre;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%