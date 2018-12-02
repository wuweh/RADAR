function [yreal,ym]=funtrackingsnake(a,omig,t,R)
yreal=a*sin(omig*t);
ym=yreal+randn(size(yreal))*sqrt(R);