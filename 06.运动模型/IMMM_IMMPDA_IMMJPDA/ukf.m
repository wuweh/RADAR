
function [X,P,e,S]=ukf(ffun,X,P,hfun,Z,Q,R)
    L=numel(X);  
    m=numel(Z);  

    alpha=1e-2;
    ki=0;
    beta=2;

    lambda=alpha^2*(L+ki)-L; 
    c=L+lambda;              
    Wm=[lambda/c 0.5/c+zeros(1,2*L)]; 
    Wc=Wm;
    Wc(1)=Wc(1)+(1-alpha^2+beta);

    c=sqrt(c);
    Xsigmaset=sigmas(X,P,c); 
    [X1means,X1,P1,X2]=ut(ffun,Xsigmaset,Wm,Wc,L,Q);   
    [Zpre,Z1,S,Z2]=ut(hfun,X1,Wm,Wc,m,R);

    Pxz=X2*diag(Wc)*Z2';
    K=Pxz*inv(S); 
    e = Z-Zpre;
    X=X1means+K*e;   
    P=P1-K*Pxz';           
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sigma points around reference point
% Inputs:
%   x: reference point
%   P: covariance
%   c: coefficient
% Output:
%   X: Sigma points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xmeans,Xsigma_pre,P,Xdiv]=ut(fun,Xsigma,Wm,Wc,n,COV)
    LL=size(Xsigma,2);
    Xmeans=zeros(n,1);
    Xsigma_pre=zeros(n,LL);
    for k=1:LL
        Xsigma_pre(:,k)=fun*Xsigma(:,k);  
        Xmeans=Xmeans+Wm(k)*Xsigma_pre(:,k); 
    end
    Xdiv=Xsigma_pre-Xmeans(:,ones(1,LL));
    P=Xdiv*diag(Wc)*Xdiv'+COV;              
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sigma points around reference point
% Inputs:
%   x: reference point
%   P: covariance
%   c: coefficient
% Output:
%   X: Sigma points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xset=sigmas(X,P,c)
    A = c*chol(P)';
    Y = X(:,ones(1,numel(X)));
    Xset = [X Y+A Y-A];
end