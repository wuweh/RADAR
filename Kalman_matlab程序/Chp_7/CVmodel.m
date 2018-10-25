function [A,Q]=CVmodel(T,qq)
A=[1 T;0 1];
Q=[T^3/3 T^2/2;T^2/2 T]*qq;