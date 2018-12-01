function [xPts, wPts, nPts] = scaledSymmetricSigmaPoints(x,P,alpha,beta,kappa)
n    = size(x(:),1);
nPts = 2*n+1;           
kappa = alpha^2*(n+kappa)-n;
wPts=zeros(1,nPts);
xPts=zeros(n,nPts);
Psqrtm=(chol((n+kappa)*P))';  
xPts=[zeros(size(P,1),1) -Psqrtm Psqrtm];
xPts = xPts + repmat(x,1,nPts);  
wPts=[kappa 0.5*ones(1,nPts-1) 0]/(n+kappa);
wPts(nPts+1) = wPts(1) + (1-alpha^2) + beta;