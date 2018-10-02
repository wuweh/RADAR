function [x,P] = ukf(fstate, x, P, hmeas, z, Q, R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UKF Unscented Kalman Filter for nonlinear dynamic systems
% ���𿨶����˲���Unscented Kalman Filter�������������ڶ�̬������ϵͳ
% for nonlinear dynamic system (noises are assumed as additive):
%   x_k+1 = f(x_k) + w_k
%   z_k = h(x_k) + v_k
% w ~ N(0,Q) meaning w is gaussian noise with covariance Q
% v ~ N(0,R) meaning v is gaussian noise with covariance R
% =============================����˵��=================================
% Inputs: 
% fstate  -[function]: ״̬����f(x)
%     x   -     [vec]: ״̬������� "a priori" state estimate
%     P   -     [mat]: ����������� "a priori" estimated state covariance
% hmeas   -[function]: ���ⷽ��h(x)
%     z   -     [vec]: ��������     current measurement
%     Q   -     [mat]: ״̬��������w(t) process noise covariance
%     R   -     [mat]: ���ⷽ������v(t) measurement noise covariance
% Output:
%     x   -     [mat]: ״̬������� "a posteriori" state estimate
%     P   -     [mat]: ���������� "a posteriori" state covariance
% =====================================================================
% By Yi Cao at Cranfield University, 04/01/2008
% Modified by JD Liu 2010-4-20
% Modified by zhangwenyu, 12/23/2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<7
    error('Not enough inputarguments!');
end

% ��ʼ����Ϊ�˼򻯺�������lamda�Ĺ��̱�Ĭ��
L = numel(x);                                 %numer of states ��6
m = numel(z);                                 %numer of measurements 6
alpha = 1e-3;                                 %default, tunable
ki = 0;                                       %default, tunable
beta = 2;                                     %default, tunable

% UTת������
lambda = alpha^2*(L+ki)-L;                    %scaling factor
c = L+lambda;                                 %scaling factor
Wm = [lambda/c 0.5/c+zeros(1,2*L)];           %weights for means Ȩ�ؼ���
Wc = Wm;
Wc(1) = Wc(1)+(1-alpha^2+beta);               %weights for covariance
c = sqrt(c);

X = sigmas(x,P,c);                            %sigma points around x �����õ�sigmaֵ
[x1,X1,P1,X2] = ut(fstate,X,Wm,Wc,L,Q);       %unscented transformation of process
% X1=sigmas(x1,P1,c);                         %sigma points around x1
% X2=X1-x1(:,ones(1,size(X1,2)));             %deviation of X1
[z1,Z1,P2,Z2] = ut(hmeas,X1,Wm,Wc,m,R);       %unscented transformation of measurments

% �˲�����
%��Э����
P12 = X2*diag(Wc)*Z2';                        %transformed cross-covariance
%�˲�����
K = P12*inv(P2);   
%����״̬����
x = x1+K*(z-z1);                              %state update
%Э�������
P = P1-K*P12';                                %covariance update

