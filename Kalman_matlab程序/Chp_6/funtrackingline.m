function [yreal,ym]=funtrackingline(a,t,R)
yreal=a*t;
ym=yreal+randn(size(yreal))*sqrt(R);