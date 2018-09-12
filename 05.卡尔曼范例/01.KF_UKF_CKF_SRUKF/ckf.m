
function [x, P] = ckf(xhat, Pplus, z)
global Q R fai gama kesi w m;
%% ---------------------------3rd-degree CKF----------------------
%% -----------------------------Time Update-----------------------
% Evaluate the Cholesky factor
Shat = chol(Pplus, 'lower');

for cpoint = 1 : m
    % Evaluate the cubature points
    rjpoint(:, cpoint) = xhat + Shat*kesi(:,cpoint);
    % Evaluate the propagated cubature points
    Xminus(:, cpoint) = fai(rjpoint(:, cpoint));
end

% Estimate the redicted state
xminus_sum = Xminus(:,1);
for k = 2:m
    xminus_sum = xminus_sum + Xminus(:,k);
end

xminus_sum = w*xminus_sum;

% Estimate the predicted error covariance
a = zeros(6,1);
for k = 1:m
    a = a+(Xminus(:,k)-xminus_sum) * ((Xminus(:,k)-xminus_sum)');
end
Pminus = w*a + Q;
%% ---------------------------------------------------------------

%% -------------------------Measurement Update--------------------
% Evaluate the Cholesky factor
Sminus = chol(Pminus, 'lower');
for cpoint = 1 : m
    % Evaluate the cubature points
    rjpoint1(:, cpoint) = xminus_sum + Sminus*kesi(:, cpoint);
    % Evaluate the propagated cubature points
%     h=@(x)[ sqrt(x(1)^2+x(2)^2);
%         atan(x(2)/x(1))];
    Z(:,cpoint) = gama(rjpoint1(:,cpoint));
end

zhat = Z(:,1);
for k = 2:m
    zhat = zhat + Z(:,k);
end

zhat_sum = w*zhat;

% Estimate the predicted error covariance
a = zeros(2,1);
b = zeros(6,2);

for k = 1:m
    a = a+(Z(:,k)-zhat_sum) * (Z(:,k)-zhat_sum)';
    b = b+(Xminus(:,k)-xminus_sum) * (Z(:,k)-zhat_sum)';
end

Pzminus = w*a + R;
Pxzminus = w*b;

% Estimate the Kalman gain
W = Pxzminus * inv(Pzminus);

% Estimate the updated state
x = xminus_sum + W * (z - zhat_sum);

% Estimate the correspondong error covariance
P = Pminus - W * Pzminus * W';
end