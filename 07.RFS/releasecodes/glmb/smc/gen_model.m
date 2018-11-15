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
model.P_S= .99;
model.Q_S= 1-model.P_S;

% birth parameters (LMB birth model, single component only)
model.T_birth= 4;         %no. of LMB birth terms
model.L_birth= zeros(model.T_birth,1);                                          %no of Gaussians in each LMB birth term
model.r_birth= zeros(model.T_birth,1);                                          %prob of birth for each LMB birth term
model.w_birth= cell(model.T_birth,1);                                           %weights of GM for each LMB birth term
model.m_birth= cell(model.T_birth,1);                                           %means of GM for each LMB birth term
model.B_birth= cell(model.T_birth,1);                                           %std of GM for each LMB birth term
model.P_birth= cell(model.T_birth,1);                                           %cov of GM for each LMB birth term

model.L_birth(1)=1;                                                             %no of Gaussians in birth term 1
model.r_birth(1)=0.02;                                                          %prob of birth
model.w_birth{1}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{1}(:,1)= [ -1500; 0; 250; 0; 0 ];                                 %mean of Gaussians
model.B_birth{1}(:,:,1)= diag([ 50; 50; 50; 50; 6*(pi/180) ]);                  %std of Gaussians
model.P_birth{1}(:,:,1)= model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';      %cov of Gaussians

model.L_birth(2)=1;                                                             %no of Gaussians in birth term 2
model.r_birth(2)=0.02;                                                          %prob of birth
model.w_birth{2}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{2}(:,1)= [ -250; 0; 1000; 0; 0 ];                                 %mean of Gaussians
model.B_birth{2}(:,:,1)= diag([ 50; 50; 50; 50; 6*(pi/180) ]);                  %std of Gaussians
model.P_birth{2}(:,:,1)= model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';      %cov of Gaussians

model.L_birth(3)=1;                                                             %no of Gaussians in birth term 3
model.r_birth(3)=0.03;                                                          %prob of birth
model.w_birth{3}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{3}(:,1)= [ 250; 0; 750; 0; 0 ];                                   %mean of Gaussians
model.B_birth{3}(:,:,1)= diag([ 50; 50; 50; 50; 6*(pi/180) ]);                  %std of Gaussians
model.P_birth{3}(:,:,1)= model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';      %cov of Gaussians

model.L_birth(4)=1;                                                             %no of Gaussians in birth term 4
model.r_birth(4)=0.03;                                                          %prob of birth
model.w_birth{4}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{4}(:,1)= [ 1000; 0; 1500; 0; 0 ];                                 %mean of Gaussians
model.B_birth{4}(:,:,1)= diag([ 50; 50; 50; 50; 6*(pi/180) ]);                  %std of Gaussians
model.P_birth{4}(:,:,1)= model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';      %cov of Gaussians

% observation model parameters (noisy r/theta only)
% measurement transformation given by gen_observation_fn, observation matrix is N/A in non-linear case
model.D= diag([ 2*(pi/180); 10 ]);      %std for angle and range noise
model.R= model.D*model.D';              %covariance for observation noise

% detection parameters
model.P_D= .98;   %probability of detection in measurements
model.Q_D= 1-model.P_D; %probability of missed detection in measurements

% clutter parameters
model.lambda_c= 10;                             %poisson average rate of uniform clutter (per scan)
model.range_c= [ -pi/2 pi/2; 0 2000 ];          %uniform clutter on r/theta
model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density




