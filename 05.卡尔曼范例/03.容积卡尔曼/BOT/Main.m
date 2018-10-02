%% ----------------------------------------------------------------
% Bearing-only tracking model for tesing the cubature Kalman filter (CKF) 
% and the embedded cubature Kalman filter (ECKF)
% 
% Version: 20140307
% 
% 
% @All rights reserved by Xin-Chun Zhang from UESTC
% You can use the code for learning, teaching and education, do not use it
% for commerical purpose!
% 
% 
% 
% 参考：
% [1] I.Arasaratnam, S.Haykin. Cubature Kalman Filters[J]. IEEE Trans. Automat.
%     Control, 2009, 54(6)
% [2] 张鑫春,等. 一种用于目标跟踪的嵌入式容积卡尔曼滤波器设计[C]. 综合电子技术教
%     育部重点实验室2012学术年会论文集，2012
% [3] 张鑫春, 等. "均方根嵌入式容积卡尔曼滤波," 控制理论与应用, 2013, 30(9).         
% [4] Xin-Chun Zhang, Yun-Long Teng. A new derivation of the cubature
%     Kalman filters [J]. Asian Journal of Control, 2013, Accept
% [5] Xin-Chun Zhang, A novel cubature Kalman filter for nonlinear state
%     estimation [C]. 52nd IEEE Conference on Decision and Control,
%     Florence, Italy, 2013
% [6] Zhang Xin-Chun, Guo Cheng-Jun, Cubature Kalman filters: Derivation
%     and extension [J]. Chinese Physics B, 2013, 22(12)
% [7] Xin-Chun Zhang. Cubature information filter using embedded and
%     high-degree cubature rules [J]. Circuits, Systems and Signal
%     processing, Accept(2013)/Online(2014)
%
% 
% 
% 
%
%                      By Xin-Chun Zhang 
%                 E-mail：irving_zhang@163.com
%       University of Electronic Science and Technology of China
%% ----------------------------------------------------------------

clear all;
close all;
clc;

h = waitbar(0, '1', 'Name', 'Please wait ...', 'WindowStyle', 'modal', ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
setappdata(h, 'canceling', 0);

% Global variables
global Q R fai gama kesi w n m w1 w2 kesi_im3 b ;

Q = 0.0005^2; % process noise
R = 0.001^2; % measurement noise
fai = [1 1 0 0; 0 1 0 0; 0 0 1 1; 0 0 0 1]; %状态转移矩阵
gama = [0.5 0; 1 0; 0 0.5; 0 1];
n = 4; %dimension of the system
m = 2 * n; % Number of cubature points

% 3rd-degree CKF/SCKF
w = 1 / m; % Weight
kesi1 = eye(n);
kesi2 = -eye(n);
kesi = [kesi1, kesi2] * (sqrt(m / 2)); % Construction of the cubature points

% 3rd-degree ECKF
delta = 1;
w1 = 1 / (2^(n + 1) * delta);
w2 = 1 - 1 / (2 * delta);
yitao = zeros(n, 1);
yitak1 = [1 1 1 1; 1 1 1 -1; 1 1 -1 -1; 1 -1 -1 -1; 1 -1 1 1; 1 -1 -1 1; ...
    1 -1 1 -1; 1 1 -1 1]' * sqrt(2 * delta);
kesi_im3 = [yitak1, -yitak1, yitao];
[a, b] = size(kesi_im3);

num = 30; % steps of the simulation

xarray_sum = zeros(n, num + 1);

xhatarray_sum = zeros(n, num + 1);
error_avg = zeros(n, num + 1);

xhatarray_S_sum = zeros(n, num + 1);
error_S_avg = zeros(n, num + 1);

xhatarray_5_sum = zeros(n, num + 1);
error_5_avg = zeros(n, num + 1);

xhatarray_I_sum = zeros(n, num + 1);
error_I_avg = zeros(n, num + 1);

len = 100;
for iter = 1 : len
    if getappdata(h, 'canceling')
        iter = iter - 1;
        break;
    end
    waitbar(iter / len, h, sprintf('%s%%', num2str(100 * iter / len)));
    x = [-0.05, 0.0001, 0.7, -0.055]'; % actual position
%     Estimates position
    xhat = x;  
    xhat_5 = x;
    
    Pplus = diag([0.01, 0.005^2, 0.01, 0.01^2]);
    Pplus_5 = Pplus;
    Splus = Pplus;
    
    xarray = [x];
    xhatarray = [x];
    xhatarray_5 = [x];
    
    for i = 1 : num
        x = fai * x + gama * sqrt(Q) * randn(2, 1);  % Process equation
        z = atan(x(3) / x(1)) + sqrt(R) * randn;  % measurement equation
        
        [xhat, Splus] = Third_degree_ECKF(xhat, Splus, z);
        [xhat_5, Pplus_5] = Third_degree_CKF(xhat_5, Pplus_5, z);
        
        xarray = [xarray x];
        xhatarray = [xhatarray xhat];
        xhatarray_5 = [xhatarray_5 xhat_5];
    end
    %% ---------------------------------------------------------------
    xarray_sum = xarray_sum + xarray;

    xhatarray_sum = xhatarray_sum + xhatarray;

    xhatarray_5_sum = xhatarray_5_sum + xhatarray_5;
    
end
delete(h);
xarray_avg = xarray_sum / iter;

xhatarray_avg = xhatarray_sum / iter;

xhatarray_5_avg = xhatarray_5_sum / iter;

k = 0 : num;

figure('Units', 'Normalized', 'Color', [1 1 1], 'Name', 'Trajectory', ...
    'NumberTitle', 'off', 'Position', [0.1 0.1 0.8 0.8], 'Tag', 'newplot'), hold on;
box on;
plot(xarray_avg(1, :), xarray_avg(3, :), 'b-x', xhatarray_avg(1, :), xhatarray_avg(3, :), 'r-*', ...
    xhatarray_5_avg(1, :), xhatarray_5_avg(3, :), 'k-p');
set(gca, 'fontname', 'Arial', 'fontsize', 18);
set(gcf, 'Color', 'White');
axis tight;
legend('True', 'ECKF', 'CKF', 'Orientation', 'vertical', 'Location', 'SouthWest');
legend('boxoff');
xlabel('x', 'fontname', 'Arial', 'fontsize', 22);
ylabel('y', 'fontname', 'Arial', 'fontsize', 22);
