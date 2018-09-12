function [x, P] = Third_degree_CKF(xhat, Pplus, z)
global Q R fai gama kesi w m;
%% ---------------------------3rd-degree CKF----------------------

%% -----------------------------Time Update-----------------------
% Evaluate the Cholesky factor
Shat = chol(Pplus, 'lower');
for cpoint = 1 : m
    % Evaluate the cubature points
    rjpoint(:, cpoint) = Shat * kesi(:, cpoint) + xhat;
    % Evaluate the propagated cubature points
    Xminus(:, cpoint) = fai * rjpoint(:, cpoint);
end

% Estimate the predicted state
xminus = w * sum(Xminus, 2);

% Estimate the predicted error covariance
Pminus = w * (Xminus * Xminus') - xminus * xminus' + gama * Q * gama';
%% ---------------------------------------------------------------

%% -------------------------Measurement Update--------------------
% Evaluate the Cholesky factor
Sminus = chol(Pminus, 'lower');
for cpoint = 1 : m
    % Evaluate the cubature points
    rjpoint1(:, cpoint) = Sminus * kesi(:, cpoint) + xminus;
    % Evaluate the propagated cubature points
    Z(cpoint) = atan(rjpoint1(3, cpoint) / rjpoint1(1, cpoint));
end
% Estimate the predicted measurement
zhat = w * sum(Z);

% Estimate the innovation covariance matrix
Pzminus = w * sum(Z * Z') - zhat^2 + R;

% Estimate the cross-covariance matrix
Pxzminus = w * rjpoint1 * Z' - xminus * zhat';

% Estimate the Kalman gain
W = Pxzminus * inv(Pzminus);

% Estimate the updated state
x = xminus + W * (z - zhat);

% Estimate the correspondong error covariance
P = Pminus - W * Pzminus * W';