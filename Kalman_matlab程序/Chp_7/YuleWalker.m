function [a,q]=YuleWalker(p)
a=[];q=[];
r0=0;
r1=0;
for k=2:length(p)
    r0=r0+(p(k)*p(k)-r0)/k;
    r1=r1+(p(k)*p(k-1)-r1)/k;
    P=r1;
H=r0;
a0=inv(H'*H)*H'*P;
q0=r0-a0*r1;
a=[a a0];
q=[q q0];
end