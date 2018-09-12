%% ---------------------------------------------------------------
% Third-degree embedded/imbedded cubature Kalman filter (ECKF/ICKF)
% Both the ECKF or ICKF and its square-root version (SECKF/SICKF) are 
% proposed by Xin-Chun Zhang from UESTC.
% 
% 
%% ---------------------------------------------------------------
function [x, P] = Third_degree_ECKF(xhat, Pplus, z)

global Q R fai gama m w1 w2 kesi_im3 b;

%% ------------------------3rd-degree ECKF------------------------

%% --------------------------Time Update--------------------------
% Evaluate the Cholesky factor
Shat = chol(Pplus, 'lower');
for cpoint = 1 : b
    % Evaluate the cubature points
    rjpoint(:, cpoint) = Shat * kesi_im3(:, cpoint) + xhat;
    % Evaluate the propagated cubature points
    Xminus(:, cpoint) = fai * rjpoint(:, cpoint);
end


% Estimate the predicted state
xminus = w1 * sum(Xminus(:, 1 : 2 * m), 2) + w2 * Xminus(:, b);

% Estiamate the predicted error covariance
Pminus = gama * Q * gama' + w1 * (Xminus(:, 1 : 2 * m) - repmat(xminus, 1 , 2 * m)) * (Xminus(:, 1 : 2 * m)...
    - repmat(xminus, 1, 2 * m))' + w2 * (Xminus(:, b) - xminus) * (Xminus(:, b) - xminus)';
%% ---------------------------------------------------------------

%% ------------------------Measurement Update---------------------
% Evaluate the Cholesky factor
Sminus = chol(Pminus, 'lower');
for cpoint = 1 : b
    % Evaluate the cubature points
    rjpoint1(:, cpoint) = Sminus * kesi_im3(:, cpoint) + xminus;
    % Evaluate the propagated cubature points
    Z(cpoint) = atan(rjpoint1(3, cpoint) / rjpoint1(1, cpoint));
end

% Estimate the predicted measurement
zhat = w1 * sum(Z(1 : 2 * m)) + w2 * Z(b);

% Estimate the innovation covariance matrix
Pzminus = R + w1 * (Z(1 : 2 * m) - repmat(zhat, 1, 2 * m)) * (Z(1 : 2 * m) ...
    - repmat(zhat, 1, 2 * m))' + w2 * (Z(b) - zhat) * (Z(b) - zhat)';

% Estimate the cross-covariance matrix
Pxzminus = w1 * (rjpoint1(:, 1 : 2 * m) - repmat(xminus, 1, 2 * m)) * (Z(1 : 2 * m) - ...
    repmat(zhat,1,2*m))' + w2 * (rjpoint1(:, b) - xminus) * (Z(b) - zhat)';

% Estimate the Kalman gain
W = Pxzminus * inv(Pzminus);

% Estimate the updated state
x = xminus + W * (z - zhat);

% Estimate the corresponding error covariance
P = Pminus - W * Pzminus * W';
%% ---------------------------------------------------------------