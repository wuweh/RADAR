function [a,q]=LSRecursive(p)
a0=0.5;
Q0=100;
q0=0;
a=[];q=[];
a=[a a0];
q=[q q0];
for k=2:length(p)
    q0=((k-1)*q0+(p(k)-a0*p(k-1))^2)/k;
    K=Q0*p(k-1)/(q0+p(k-1)*Q0*p(k-1));
    a0=a0+K*(p(k)-a0*p(k-1));
    Q0=(1-K*p(k-1))*Q0;
    a=[a a0];
    q=[q q0];
end