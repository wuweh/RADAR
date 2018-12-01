%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ×´Ì¬·½³Ìº¯Êý
function [y] = ffun(x,t);
 
if nargin < 2
    error('Not enough input arguments.'); 
end
beta = 0.5;
y = 1 + sin(4e-2*pi*t) + beta*x;