function [xkk1,Skk1] = Predict(xkk,Skk)

global Qsqrt;

xkk1 = xkk;

[foo,Skk1] = qr([Skk Qsqrt]',0);

Skk1 = Skk1';

