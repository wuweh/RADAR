function model= gen_model

% basic parameters
model.x_dim= 5+1;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector
model.v_dim= 3+1;   %dimension of process noise
model.w_dim= 2;   %dimension of observation noise

% dynamical model parameters (CT model)
% state transformation given by gen_newstate_fn, transition matrix is N/A in non-linear case
model.T= 1;                         %sampling period
model.sigma_vel= 5;
model.sigma_turn= (pi/180);   %std. of turn rate variation (rad/s)
model.bt= model.sigma_vel*[ (model.T^2)/2; model.T ];
model.B2= [ model.bt zeros(2,2); zeros(2,1) model.bt zeros(2,1); zeros(1,2) model.T*model.sigma_turn ];
model.B= eye(model.v_dim-1);
model.Q= model.B*model.B';
model.pdvarfac_tg= 0.01^2;

% survival/death parameters
% use compute_pS for state dependent parameterization
% use compute_qS for state dependent parameterization

% birth parameters (Poisson birth model, multiple Gaussian components)
model.L_birth= 4;                                                     %no. of Gaussian birth terms
model.w_birth= zeros(model.L_birth,1);                                %weights of Gaussian birth terms (per scan) [sum gives average rate of target birth per scan]
model.m_birth= zeros(model.x_dim-1,model.L_birth);                      %means of Gaussian birth terms 
model.B_birth= zeros(model.x_dim-1,model.x_dim-1,model.L_birth);          %std of Gaussian birth terms
model.P_birth= zeros(model.x_dim-1,model.x_dim-1,model.L_birth);          %cov of Gaussian birth terms
model.u_b = zeros(model.L_birth,1); model.v_b = zeros(model.L_birth,1);     % Beta parameters for unknown detection profile


model.w_birth(1)= 2/100;                                              %birth term 1
model.m_birth(:,1)= [ -1500; 0; 250; 0; 0 ];
model.B_birth(:,:,1)= diag([ 50; 50; 50; 50; 6*(pi/180) ]);
model.P_birth(:,:,1)= model.B_birth(:,:,1)*model.B_birth(:,:,1)';
model.u_b(1) = 95;
model.v_b(1) = 5;

model.w_birth(2)= 2/100;                                              %birth term 2
model.m_birth(:,2)= [ -250; 0; 1000; 0; 0 ];
model.B_birth(:,:,2)= diag([ 50; 50; 50; 50; 6*(pi/180) ]);
model.P_birth(:,:,2)= model.B_birth(:,:,2)*model.B_birth(:,:,2)';
model.u_b(2) = 95;
model.v_b(2) = 5;

model.w_birth(3)= 3/100;                                              %birth term 3
model.m_birth(:,3)= [ 250; 0; 750; 0; 0 ]; 
model.B_birth(:,:,3)= diag([ 50; 50; 50; 50; 6*(pi/180) ]);
model.P_birth(:,:,3)= model.B_birth(:,:,3)*model.B_birth(:,:,3)';
model.u_b(3) = 95;
model.v_b(3) = 5;

model.w_birth(4)= 3/100;                                              %birth term 4
model.m_birth(:,4)= [ 1000; 0; 1500; 0; 0 ];
model.B_birth(:,:,4)= diag([ 50; 50; 50; 50; 6*(pi/180) ]);
model.P_birth(:,:,4)= model.B_birth(:,:,4)*model.B_birth(:,:,4)';
model.u_b(4) = 95;
model.v_b(4) = 5;

% observation model parameters (noisy r/theta only)
% measurement transformation given by gen_observation_fn, observation matrix is N/A in non-linear case
model.D= diag([ 2*(pi/180); 10 ]);      %std for angle and range noise
model.R= model.D*model.D';              %covariance for observation noise

%===here we set up the state spc eqn. x= Ax_old + Bv with RW model
model.Tc= 1;   %sampling period
model.Q_c= diag([1000 0 500 0 0].^2);
model.pdvarfac_c= 0.07^2;

%=== parameters for the observation 
model.C_c= [ 0 1 0 0 0 0; 0 0 0 1 0 0 ];  
model.D_c= diag([ 20*(pi/180); 400 ]);   %std for angle and range noise
model.R_c= model.D_c*model.D_c';

% detection parameters
% use compute_pD for state dependent parameterization
% use compute_qD for state dependent parameterization

% clutter parameters
model.range_c= [ -pi/2 pi/2; 0 2000 ];          %uniform clutter on r/theta
model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density

% clutter model for unknown clutter rate
model.clutter_P_S = 0.9;    % survival probability for clutter targets
model.clutter_P_D = 0.5;     % detection probability for clutter targets
model.clutter_Nt = 20;     % number of clutter generators
model.lambda_c = model.clutter_Nt * model.clutter_P_D;

% birth parameters (Poisson birth model, multiple Gaussian components)
model.Lc_birth= 1;                                                     %no. of Gaussian birth terms
model.wc_birth= zeros(model.Lc_birth,1);                                %weights of Gaussian birth terms (per scan) [sum gives average rate of target birth per scan]
model.mc_birth= zeros(model.x_dim-1,model.Lc_birth);                      %means of Gaussian birth terms 
model.Bc_birth= zeros(model.x_dim-1,model.x_dim-1,model.Lc_birth);          %std of Gaussian birth terms
model.Pc_birth= zeros(model.x_dim-1,model.x_dim-1,model.Lc_birth);          %cov of Gaussian birth terms
model.uc_b = zeros(model.Lc_birth,1); model.vc_b = zeros(model.Lc_birth,1);     % Beta parameters for unknown detection profile

for i=1:model.Lc_birth
    model.wc_birth(i)= 1;                                              %birth term 1
    model.mc_birth(:,i)= [ 0; 0; 0; 0; 0 ];
    model.Bc_birth(:,:,i)= diag([ 2000; 0; 2000; 0; 0 ]);
    model.Pc_birth(:,:,i)= model.Bc_birth(:,:,i)*model.Bc_birth(:,:,i)';
    model.uc_b(i) = 5;
    model.vc_b(i) = 5;
end

