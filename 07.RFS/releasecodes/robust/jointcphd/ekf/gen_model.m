function model= gen_model

% basic parameters
model.x_dim= 5;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector
model.v_dim= 3;   %dimension of process noise
model.w_dim= 2;   %dimension of observation noise

% dynamical model parameters (CT model)
% state transformation given by gen_newstate_fn, transition matrix is N/A in non-linear case
model.T= 1;                         %sampling period
model.sigma_vel= 15;
model.sigma_turn= (pi/180);   %std. of turn rate variation (rad/s)
model.bt= model.sigma_vel*[ (model.T^2)/2; model.T ];
model.B2= [ model.bt zeros(2,2); zeros(2,1) model.bt zeros(2,1); zeros(1,2) model.T*model.sigma_turn ];
model.B= eye(model.v_dim);
model.Q= model.B*model.B';

% survival/death parameters
model.P_S= .99;
model.Q_S= 1-model.P_S;

% birth parameters (Poisson birth model, multiple Gaussian components)
model.lambda_b = 0.1;
model.L_birth= 4;                                                     %no. of Gaussian birth terms
model.w_birth= zeros(model.L_birth,1);                                %weights of Gaussian birth terms (per scan) [sum gives average rate of target birth per scan]
model.m_birth= zeros(model.x_dim,model.L_birth);                      %means of Gaussian birth terms 
model.B_birth= zeros(model.x_dim,model.x_dim,model.L_birth);          %std of Gaussian birth terms
model.P_birth= zeros(model.x_dim,model.x_dim,model.L_birth);          %cov of Gaussian birth terms
model.u_b = zeros(model.L_birth,1); model.v_b = zeros(model.L_birth,1);     % Beta parameters for unknown detection profile


model.w_birth(1)= 2/100;                                              %birth term 1
model.m_birth(:,1)= [ -1500; 0; 250; 0; 0 ];
model.B_birth(:,:,1)= diag([ 70; 70; 70; 70; 7*(pi/180) ]);
model.P_birth(:,:,1)= model.B_birth(:,:,1)*model.B_birth(:,:,1)';
model.u_b(1) = 98;
model.v_b(1) = 2;

model.w_birth(2)= 2/100;                                              %birth term 2
model.m_birth(:,2)= [ -250; 0; 1000; 0; 0 ];
model.B_birth(:,:,2)= diag([ 70; 70; 70; 70; 7*(pi/180) ]);
model.P_birth(:,:,2)= model.B_birth(:,:,2)*model.B_birth(:,:,2)';
model.u_b(2) = 98;
model.v_b(2) = 2;

model.w_birth(3)= 3/100;                                              %birth term 3
model.m_birth(:,3)= [ 250; 0; 750; 0; 0 ]; 
model.B_birth(:,:,3)= diag([ 70; 70; 70; 70; 7*(pi/180) ]);
model.P_birth(:,:,3)= model.B_birth(:,:,3)*model.B_birth(:,:,3)';
model.u_b(3) = 98;
model.v_b(3) = 2;

model.w_birth(4)= 3/100;                                              %birth term 4
model.m_birth(:,4)= [ 1000; 0; 1500; 0; 0 ];
model.B_birth(:,:,4)= diag([ 70; 70; 70; 70; 7*(pi/180) ]);
model.P_birth(:,:,4)= model.B_birth(:,:,4)*model.B_birth(:,:,4)';
model.u_b(4) = 98;
model.v_b(4) = 2;

% observation model parameters (noisy r/theta only)
% measurement transformation given by gen_observation_fn, observation matrix is N/A in non-linear case
model.D= diag([ 1*(pi/180); 5 ]);      %std for angle and range noise
model.R= model.D*model.D';              %covariance for observation noise

% detection parameters
model.P_D= .98;   %probability of detection in measurements
model.Q_D= 1-model.P_D; %probability of missed detection in measurements

% clutter parameters
model.range_c= [ -pi/2 pi/2; 0 2000 ];          %uniform clutter on r/theta
model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density

model.lambda_cb= 1;
model.w_cb = model.lambda_cb;
model.clutter_P_S= 0.9;
model.clutter_P_D = 0.5;     % detection probability for clutter targets
model.Lc_birth = 1; model.u_cb = 1; model.v_cb = 1;
model.clutter_Nt = 20;     % number of clutter generators
model.lambda_c = model.clutter_Nt * model.clutter_P_D;



