function [x,P] = ukf(fstate, x, P, hmeas, z, Q, R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UKF Unscented Kalman Filter for nonlinear dynamic systems
% 无损卡尔曼滤波（Unscented Kalman Filter）函数，适用于动态非线性系统
% for nonlinear dynamic system (noises are assumed as additive):
%   x_k+1 = f(x_k) + w_k
%   z_k = h(x_k) + v_k
% w ~ N(0,Q) meaning w is gaussian noise with covariance Q
% v ~ N(0,R) meaning v is gaussian noise with covariance R
% =============================参数说明=================================
% Inputs: 
% fstate  -[function]: 状态方程f(x)
%     x   -     [vec]: 状态先验估计 "a priori" state estimate
%     P   -     [mat]: 方差先验估计 "a priori" estimated state covariance
% hmeas   -[function]: 量测方程h(x)
%     z   -     [vec]: 量测数据     current measurement
%     Q   -     [mat]: 状态方程噪声w(t) process noise covariance
%     R   -     [mat]: 量测方程噪声v(t) measurement noise covariance
% Output:
%     x   -     [mat]: 状态后验估计 "a posteriori" state estimate
%     P   -     [mat]: 方差后验估计 "a posteriori" state covariance
% =====================================================================
% By Yi Cao at Cranfield University, 04/01/2008
% Modified by JD Liu 2010-4-20
% Modified by zhangwenyu, 12/23/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<7
    error('Not enough inputarguments!');
end

% 初始化，为了简化函数，求lamda的过程被默认
L = numel(x);                                 %numer of states 。6
m = numel(z);                                 %numer of measurements 6
alpha = 1e-3;                                 %default, tunable
ki = 0;                                       %default, tunable
beta = 2;                                     %default, tunable

% UT转换部分
lambda = alpha^2*(L+ki)-L;                    %scaling factor
c = L+lambda;                                 %scaling factor
Wm = [lambda/c 0.5/c+zeros(1,2*L)];           %weights for means 权重计算
Wc = Wm;
Wc(1) = Wc(1)+(1-alpha^2+beta);               %weights for covariance
c = sqrt(c);

X = sigmas(x,P,c);                            %sigma points around x 真正得到sigma值
[x1,X1,P1,X2] = ut(fstate,X,Wm,Wc,L,Q);       %unscented transformation of process
% X1=sigmas(x1,P1,c);                         %sigma points around x1
% X2=X1-x1(:,ones(1,size(X1,2)));             %deviation of X1
[z1,Z1,P2,Z2] = ut(hmeas,X1,Wm,Wc,m,R);       %unscented transformation of measurments

% 滤波部分
%互协方差
P12 = X2*diag(Wc)*Z2';                        %transformed cross-covariance
%滤波增益
K = P12*inv(P2);   
%后验状态估计
x = x1+K*(z-z1);                              %state update
%协方差矩阵
P = P1-K*P12';                                %covariance update

