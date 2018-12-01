%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 观测方程函数
function [y] = hfun(x,t);
 
if nargin < 2, 
    error('Not enough input arguments.');
end
if t<=30
    y = (x.^(2))/5;
else
    y = -2 + x/2;
end;