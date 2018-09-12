function [xkk,Skk] = Update2(xkk1,Skk1,z)
global Rsqrt;
%%%========================================================================
%%% Althought algebraically equivalent to Update.m, Update2 is
%%% mathematically more elegant (Read `Hybrid CKF' downloadble from 
%%% http://grads.ece.mcmaster.ca/~aienkaran/)
%%%========================================================================
%%%========================================================================
%%% Genrate a set of Cubature Points
%%%========================================================================

nx = 2; %state vector dimension
nz = 2; %mst vector dimension
nPts = 2*nx;
CPtArray = sqrt(nPts/2)*[eye(nx) -eye(nx)];
%%%========================================================================

Xi =  repmat(xkk1,1,nPts) + Skk1*CPtArray;
    
Zi = MstEq(Xi);
    
zkk1 = sum(Zi,2)/nPts;   % predicted Measurement

X = (Xi-repmat(xkk1,1,nPts))/sqrt(nPts);
    
Z = (Zi-repmat(zkk1,1,nPts))/sqrt(nPts);  

[foo,S] = qr([Z Rsqrt; X zeros(nx,nz)]',0);

S = S';

A = S(1:nz,1:nz);   % Square-root Innovations Covariance

B = S(nz+1:end,1:nz);

C = S(nz+1:end,nz+1:end);

G = B/A;          % Cubature Kalman Gain  

xkk = xkk1 + G*(z-zkk1);  
    
Skk = C;


