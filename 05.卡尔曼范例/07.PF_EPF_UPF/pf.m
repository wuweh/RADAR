% 基本粒子滤波算法
function [Xo,Xoset]=pf(Xiset,Z,N,k,R,g1,g2)
 
tic
 
resamplingScheme=1;
 
Zpre=ones(1,N);   
Xsetpre=ones(1,N);  
w = ones(1,N);     

 
for i=1:N
    Xsetpre(i) = feval('ffun',Xiset(i),k) + gengamma(g1,g2);
end

 
for i=1:N
    Zpre(i) = feval('hfun',Xsetpre(i),k);
    w(i) = inv(sqrt(R)) * exp(-0.5*inv(R)*((Z-Zpre(i))^(2))) ...
        + 1e-99; 
end
w = w./sum(w);             
 
if resamplingScheme == 1
    outIndex = residualR(1:N,w');
elseif resamplingScheme == 2
    outIndex = systematicR(1:N,w');  
else
    outIndex = multinomialR(1:N,w');
end

Xoset = Xsetpre(outIndex); 
Xo=mean(Xoset);