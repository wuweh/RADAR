function model= gen_model

% basic parameters
model.x_dim= 5;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector
model.v_dim= 3;   %dimension of process noise
model.w_dim= 2;   %dimension of observation noise

% dynamical model parameters (CT model)
% state transformation given by gen_newstate_fn, transition matrix is N/A in non-linear case
model.T= 1;                         %sampling period
model.sigma_vel= 5;
model.sigma_turn= (pi/180);   %std. of turn rate variation (rad/s)
model.bt= model.sigma_vel*[ (model.T^2)/2; model.T ];
model.B2= [ model.bt zeros(2,2); zeros(2,1) model.bt zeros(2,1); zeros(1,2) model.T*model.sigma_turn ];
model.B= eye(model.v_dim);
model.Q= model.B*model.B';

% survival/death parameters
% N/A for single target

% birth parameters
% N/A for single target

% observation model parameters (noisy r/theta only)
% measurement transformation given by gen_observation_fn, observation matrix is N/A in non-linear case
model.D= diag([ 2*(pi/180); 10 ]);      %std for angle and range noise
model.R= model.D*model.D';              %covariance for observation noise

% detection parameters
% use compute_pD for state dependent parameterization
% use compute_qD for state dependent parameterization

% clutter parameters
model.lambda_c= 20;                             %poisson average rate of uniform clutter (per scan)
model.range_c= [ -pi/2 pi/2; 0 2000 ];          %uniform clutter on r/theta
model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density




