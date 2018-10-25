function [A,Q]=CAmodel(T,qq)
A=[1 T T^2/2;0 1 T;0 0 1];
Q=[T^5/20 T^4/8 T^3/6; T^4/8 T^3/3 T^2/2;T^3/6 T^2/2 T]*qq;