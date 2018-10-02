function [o,h,g,q,perr,xk]=lmodeinitial(T,r,zx,zy,vxks,vyks,perr2)
o=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];
h=[1,0,0,0;0,0,1,0];
g=[T^2/2,0;T,0;0,T^2/2;0,T];
q=[10000,0;0,10000];
xk=[zx(r-1);vxks(r-1);zy(r-1);vyks(r-1)];
perr=perr2(1:4,1:4);
    
       