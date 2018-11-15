function model= gen_model

% basic parameters
model.x_dim= 5+1;   %dimension of state vector
model.z_dim= 2;   %dimension of observation vector
model.v_dim= 3+1;   %dimension of process noise
model.w_dim= 2;   %dimension of observation noise

% dynamical model parameters (CT model)
% state transformation given by gen_newstate_fn, transition matrix is N/A in non-linear case
model.T= 1;                         %sampling period
model.sigma_vel= 15;
model.sigma_turn= 3*(pi/180);   %std. of turn rate variation (rad/s)
model.bt= model.sigma_vel*[ (model.T^2)/2; model.T ];
model.B2= [ model.bt zeros(2,2); zeros(2,1) model.bt zeros(2,1); zeros(1,2) model.T*model.sigma_turn ];
model.B= eye(model.v_dim-1);
model.Q= model.B*model.B';
model.pdvarfac_tg= 0.01^2;

% observation model parameters (noisy r/theta only)
% measurement transformation given by gen_observation_fn, observation matrix is N/A in non-linear case
model.D= diag([ 1*(pi/180); 5 ]);      %std for angle and range noise
model.R= model.D*model.D';              %covariance for observation noise

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
model.u_b_tg = zeros(model.T_birth,1); 
model.v_b_tg= zeros(model.T_birth,1);

model.b_sf= 0.7;
model.L_birth(1)=1;                                                             %no of Gaussians in birth term 1
model.r_birth(1)=0.02;                                                          %prob of birth
model.w_birth{1}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{1}(:,1)= [ -1500; 0; 250; 0; 0 ];                                 %mean of Gaussians
model.B_birth{1}(:,:,1)= diag(model.b_sf*[ 100; 100; 100; 100; 10*(pi/180) ]);                  %std of Gaussians
model.P_birth{1}(:,:,1)= model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';      %cov of Gaussians
model.u_b_tg(1)= 95;
model.v_b_tg(1)= 5;

model.L_birth(2)=1;                                                             %no of Gaussians in birth term 2
model.r_birth(2)=0.02;                                                          %prob of birth
model.w_birth{2}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{2}(:,1)= [ -250; 0; 1000; 0; 0 ];                                 %mean of Gaussians
model.B_birth{2}(:,:,1)= diag(model.b_sf*[ 100; 100; 100; 100; 10*(pi/180) ]);                 %std of Gaussians
model.P_birth{2}(:,:,1)= model.B_birth{2}(:,:,1)*model.B_birth{2}(:,:,1)';      %cov of Gaussians
model.u_b_tg(2)= 95;
model.v_b_tg(2)= 5;

model.L_birth(3)=1;                                                             %no of Gaussians in birth term 3
model.r_birth(3)=0.03;                                                          %prob of birth
model.w_birth{3}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{3}(:,1)= [ 250; 0; 750; 0; 0 ];                                   %mean of Gaussians
model.B_birth{3}(:,:,1)= diag(model.b_sf*[ 100; 100; 100; 100; 10*(pi/180) ]);                  %std of Gaussians
model.P_birth{3}(:,:,1)= model.B_birth{3}(:,:,1)*model.B_birth{3}(:,:,1)';      %cov of Gaussians
model.u_b_tg(3)= 95;
model.v_b_tg(3)= 5;

model.L_birth(4)=1;                                                             %no of Gaussians in birth term 4
model.r_birth(4)=0.03;                                                          %prob of birth
model.w_birth{4}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{4}(:,1)= [ 1000; 0; 1500; 0; 0 ];                                 %mean of Gaussians
model.B_birth{4}(:,:,1)= diag(model.b_sf*[ 100; 100; 100; 100; 10*(pi/180) ]);                  %std of Gaussians
model.P_birth{4}(:,:,1)= model.B_birth{4}(:,:,1)*model.B_birth{4}(:,:,1)';      %cov of Gaussians
model.u_b_tg(4)= 95;
model.v_b_tg(4)= 5;


% detection parameters
model.P_D= .98;   %probability of detection in measurements
model.Q_D= 1-model.P_D; %probability of missed detection in measurements
model.clutter_P_D= .5;
model.clutter_N_T= 20;
model.lambda_c= model.clutter_N_T.*model.clutter_P_D;
model.range_c= [ -pi/2 pi/2; 0 2000 ];          %uniform clutter on r/theta
model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density


model.lambda_c= 10;                             %poisson average rate of uniform clutter (per scan)

model.P_S_clt= .90;
model.P_D_clt= .50;

model.T_clt= 1;
model.Q_clt= diag([1000 0 500 0 0].^2);
model.pdvarfac_clt= 0.07^2;

% clutter parameters
model.D_clt = diag([20*(pi/180); 400]); %std for angle and range noise
model.R_clt = model.D_clt*model.D_clt';


model.T_birth_clt= 20;
model.L_birth_clt= zeros(model.T_birth_clt,1);
model.r_birth_clt= zeros(model.T_birth_clt,1);                                          %prob of birth for each LMB birth term
model.lambda_b_clt= cell(model.T_birth_clt,1);
model.m_birth_clt= cell(model.T_birth_clt,1);                                           %means of GM for each LMB birth term
model.B_birth_clt= cell(model.T_birth_clt,1);                                           %std of GM for each LMB birth term
model.P_birth_clt= cell(model.T_birth_clt,1);                                           %cov of GM for each LMB birth term
model.u_b_clt= zeros(model.T_birth_clt,1); model.v_b_clt= zeros(model.T_birth_clt,1);

for i=1:model.T_birth_clt
    model.L_birth_clt(i)= 1;
    model.r_birth_clt(i)= 0.1;
    model.lambda_b_clt{i}(1)=1;
    model.m_birth_clt{i}(:,1)= [0; 0; 0; 0; 0];
    model.B_birth_clt{i}(:,:,1)= diag([2000; 0; 2000; 0; 0;]);
    model.u_b_clt(i)= 5;
    model.v_b_clt(i)= 5;
    
    for g=1:model.L_birth_clt(i)
        model.P_birth_clt{i}(:,:,g)= model.B_birth_clt{i}(:,:,g)*model.B_birth_clt{i}(:,:,g)';
    end
end
    
    
    
