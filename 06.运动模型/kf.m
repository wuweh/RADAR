function [x,p] = kf(ffun,X,P,hfun,Z,Q,R)
%KF 此处显示有关此函数的摘要
%   此处显示详细说明
    sizeofdim = size(X,1);          %计算维度
    x1 = ffun*X;                    %计算一步预测估计值
    P1 = ffun*P*ffun'+Q;            %计算预测滤波协方差矩阵
    s = hfun*P1*hfun'+R;
    gain = P1*hfun'/s;  %计算增益矩阵
    
    x = x1+gain*(Z-hfun*x1);            %计算估计值
    p =(eye(sizeofdim)-gain*hfun)*P1;   %算协方差更新方
end

