function [y,Y,P,Y1] = ut(f,X,Wm,Wc,n,R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unscented Transformation
% UT转换函数
% Input:
% f: nonlinear map
% X: sigma points
% Wm: weights for mean
% Wc: weights for covraiance
% n: numer of outputs of f
% R: additive covariance
% Output:
% y: transformed mean
% Y: transformed smapling points
% P: transformed covariance
% Y1: transformed deviations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = size(X,2);  %返回X的列数
y = zeros(n,1);
Y = zeros(n,L);
for k=1:L
    Y(:,k) = f(X(:,k));
    y = y+Wm(k)*Y(:,k);
end
Y1 = Y-y(:,ones(1,L)); %差值
P = Y1*diag(Wc)*Y1'+R;