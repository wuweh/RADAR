function X = sigmas(x,P,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sigma points around reference point
% 构造2n+1个sigma点
% Inputs:
% x: reference point
% P: covariance
% c: coefficient
% Output:
% X: Sigma points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = c*chol(P)';
Y = x(:,ones(1,numel(x)));
X = [x Y+A Y-A];