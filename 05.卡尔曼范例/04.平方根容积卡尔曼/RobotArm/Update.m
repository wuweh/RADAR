function [xkk,Skk] = Update(xkk1,Skk1,z)

global Rsqrt;

%%%========================================================================
%%% Genrate a set of Cubature Points
%%%========================================================================

nx = 2; % state vector dimension

nPts = 2*nx;        % No. of Cubature Points

CPtArray = sqrt(nPts/2)*[eye(nx) -eye(nx)];

%%%========================================================================

Xi =  repmat(xkk1,1,nPts) + Skk1*CPtArray;
    
Zi = MstEq(Xi);
    
zkk1 = sum(Zi,2)/nPts;      % predicted Measurement

X = (Xi-repmat(xkk1,1,nPts))/sqrt(nPts);
    
Z = (Zi-repmat(zkk1,1,nPts))/sqrt(nPts);  
   
[foo,Szz] = qr([Z Rsqrt]',0);  

Szz = Szz';                 % Square-root Innovations Covariance

Pxz = X*Z';
    
G = (Pxz/Szz')/Szz;         % Cubature Kalman Gain  
    
xkk = xkk1 + G*(z - zkk1);  
    
[foo,Skk] = qr([(X - G*Z)  G*Rsqrt]',0);

Skk = Skk'; 


