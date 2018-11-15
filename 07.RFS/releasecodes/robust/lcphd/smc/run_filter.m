function est = run_filter(model,meas)

% This a the MATLAB code for the SMC implementation of the lambda-CPHD filter proposed in
% R. Mahler, B.-T. Vo and B.-N. Vo, "CPHD Filtering with unknown clutter rate and detection profile," IEEE Transactions on Signal Processing, Vol. 59, No. 8, pp. 3497-3513, 2011.
% http://ba-ngu.vo-au.com/vo/MVV_PHDrobust.pdf

% There are three versions of filters for this paper
% 1) Lambda-CPHD, 2) pD-CPHD, and 3) Lambda-pD-CPHD as sequentially described in the paper.
% This is the code for Lambda-CPHD

% ---BibTeX entry
% @ARTICLE{RobustCPHD,
% author={R. P. S. Mahler and B.-T. Vo and B.-N. Vo},
% journal={IEEE Transactions on Signal Processing},
% title={CPHD Filtering With Unknown Clutter Rate and Detection Profile},
% year={2011},
% month={August},
% volume={59},
% number={8},
% pages={3497-3513}} 
%---
% 
% based on the SMC implementation of the PHD filter given in 
%
% B.-N. Vo, S. Singh and A. Doucet, "Sequential Monte Carlo methods for Bayesian Multi-target filtering with Random Finite Sets," IEEE Trans. Aerospace and Electronic Systems, Vol. 41, No. 4, pp. 1224-1245, 2005.
% http://ba-ngu.vo-au.com/vo/VSD_SMCRFS_AES05.pdf
% ---BibTeX entry
% @ARTICLE{SMCRFS,
% author={B.-N. Vo and S. Singh and A. Doucet},
% journal={IEEE Transactions on Aerospace and Electronic Systems},
% title={Sequential Monte Carlo methods for multitarget filtering with random finite sets},
% year={2005},
% month={Oct},
% volume={41},
% number={4},
% pages={1224-1245}} 
%---


%=== Setup

%output variables
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);
est.L= cell(meas.K,1);

%filter parameters
filter.J_max= 100000;                                         %total number of particles
filter.J_target= 3000;                                        %generated number of particles per expected target
filter.J_birth= model.L_birth*filter.J_target;                %generated number of particles from birth intensity

filter.N_max= 300;                   %maximum cardinality number (for cardinality distribution)

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
w_update= eps;
m_init= [0.1;0;0.1;0;0.01];
P_init= diag([100 10 100 10 1]).^2;
x_update= m_init;
Nc_init = round((size(meas.Z{1},2)-compute_pD(model,model.m_birth)*model.w_birth)/model.clutter_P_D);
cdn_update = [zeros(Nc_init-1,1);1;zeros(filter.N_max-Nc_init+1,1)];
Nc_update= Nc_init;

%recursive filtering
for k=1:meas.K
    %---intensity prediction 
    pS_vals= compute_pS(model,x_update); pS_vals= pS_vals(:);
    pS_vals_clt= model.clutter_P_S*ones(size(x_update,2),1);

    x_predict = gen_newstate_fn(model,x_update,'noise');                                                                %surviving particles
    w_predict = pS_vals.*w_update;                                                                                      %surviving weights

    x_predict= cat(2,gen_gms(model.w_birth,model.m_birth,model.P_birth,filter.J_birth),x_predict);                      %append birth particles
    w_predict= cat(1,sum(model.w_birth)*ones(filter.J_birth,1)/filter.J_birth,w_predict);                               %append birth weights
    
    %---cardinality prediction 
    %surviving cardinality distribution
    survive_cdn_predict = zeros(filter.N_max+1,1);
    survival_factor= sum(w_update.*pS_vals)/(sum(w_update)+Nc_update) + Nc_update/(sum(w_update)+Nc_update)*model.clutter_P_S;
    if isnan(survival_factor), survival_factor=0; end %catch the degernerate zero case
    q_survival_factor= 1-survival_factor;
    for j=0:filter.N_max
        idxj=j+1;
        terms= zeros(filter.N_max+1,1);
        for ell=j:filter.N_max
            idxl= ell+1;
            terms(idxl) = exp(sum(log(1:ell))-sum(log(1:j))-sum(log(1:ell-j))+j*log(survival_factor)+(ell-j)*log(q_survival_factor))*cdn_update(idxl);
        end
        survive_cdn_predict(idxj) = sum(terms);
    end
    
    %predicted cardinality= convolution of birth and surviving cardinality distribution
    cdn_predict = zeros(filter.N_max+1,1);
    lambda_b = sum(model.w_birth)+model.lambda_cb;
    for n=0:filter.N_max
        idxn=n+1;
        terms= zeros(filter.N_max+1,1);
        for j=0:n
            idxj= j+1;
            terms(idxj)= exp(-sum(lambda_b)+(n-j)*log(lambda_b)-sum(log(1:n-j)))*survive_cdn_predict(idxj);
        end
        cdn_predict(idxn) = sum(terms);
    end
    Nc_predict = model.lambda_cb + model.clutter_P_S * Nc_update;

    %normalize predicted cardinality distribution
    cdn_predict = cdn_predict/sum(cdn_predict);
        
    %---intensity update
    %number of measurements
    m= size(meas.Z{k},2);
    pD_vals= compute_pD(model,x_predict); pD_vals= pD_vals(:);
    qD_vals= 1-pD_vals;
    
    %pre calculation for likelihood values
    if m~=0
        meas_likelihood= zeros(length(w_predict),m);
        for ell=1:m
            meas_likelihood(:,ell)= compute_likelihood(model,meas.Z{k}(:,ell),x_predict)';
        end
    end
    %There is no elementary symmetric functions in lambda-CPHD filter
    missed_factor= sum(w_predict.*qD_vals)/(sum(w_predict)+Nc_predict) + Nc_predict/(sum(w_predict)+Nc_predict)*(1-model.clutter_P_D);
    if isnan(missed_factor), missed_factor=0; end %catch the degernerate zero case
    terms0 = zeros(filter.N_max+1,1);
    for n=0:filter.N_max
        idxn= n+1;
        if n < m
            terms0(idxn) = eps(0);
        else
            terms0(idxn) =  exp(sum(log(1:n))-sum(log(1:n-m))+(n-m)*log(missed_factor));
        end   
    end
    Upsilon0 = terms0;
    
    terms1 = zeros(filter.N_max,1);
    for n=0:filter.N_max
        idxn= n+1;
        if n < m+1
            terms1(idxn) = eps(0);
        else
            terms1(idxn) =  exp(sum(log(1:n))-sum(log(1:n-(m+1)))+(n-(m+1))*log(missed_factor));
        end   
    end
    Upsilon1 = terms1;
    
    %missed detection weight
    
    w_update = 1/(sum(w_predict)+Nc_predict)*(qD_vals)*(Upsilon1'*cdn_predict)/(Upsilon0'*cdn_predict).*w_predict;
    x_update = x_predict;
    
    % detection terms
    clut_w_term= zeros(m,1); w_cumsum= sum(w_update);
    for ell=1:m        
       %calculate weights first 
       K_t= Nc_predict*model.clutter_P_D/prod(model.range_c(:,2)-model.range_c(:,1));
       w_t = pD_vals.*w_predict.*meas_likelihood(:,ell);
       C= K_t + sum(w_t);
       w_t= w_t/C; w_cumsum= w_cumsum + sum(w_t);
       clut_w_term(ell)= K_t/C;
       w_update= [w_update; w_t];
       x_update= [x_update x_predict];         
    end
    Nc_update= Nc_predict/(sum(w_predict)+Nc_predict)*(1-model.clutter_P_D)*(Upsilon1'*cdn_predict)/(Upsilon0'*cdn_predict) + sum(clut_w_term);
    w_update= w_cumsum*w_update/sum(w_update);    

    
    % cardinality update
    for n=0:filter.N_max
        idxn= n+1;
        if n < m
            cdn_update(idxn) = 0;
        else
            cdn_update(idxn) =  cdn_predict(idxn)*Upsilon0(idxn)/(Upsilon0'*cdn_predict);
        end   
    end
    cdn_update = cdn_update/sum(cdn_update);
    
    
    %---for diagnostics
    w_posterior= w_update;
    
    %---resampling
    J_rsp= min(ceil(sum(w_update)*filter.J_target),filter.J_max);
    idx= randsample(length(w_update),J_rsp,true,w_update); %idx= resample(w_update,J_rsp);
    w_update= sum(w_update)*ones(J_rsp,1)/J_rsp;
    x_update= x_update(:,idx);
 
    %--- state extraction   
    if sum(w_update) > .5
        [x_c,I_c]= our_kmeans(x_update,w_update,1);
        est.N(k)= 0;
        for j=1:size(x_c,2);
            if sum(w_update(I_c{j})) > .5,
                est.X{k}= [ est.X{k} x_c(:,j) ];
                est.N(k)= est.N(k)+1;
                est.L{k}= [];
            end
        end
    else
        est.N(k)= 0; est.X{k}= [];
    end
    est.L{k} = Nc_update*model.clutter_P_D;
    
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #eap lambda=' num2str(est.L{k},4),...
         ' #eap target=' num2str(sum(w_update),4),...
         ' #est card=' num2str(est.N(k),4),...
         ' Neff_updt= ',num2str(round(1/sum((w_posterior/sum(w_posterior)).^2)))...
         ' Neff_rsmp= ',num2str(round(1/sum((w_update/sum(w_update)).^2)))   ]);
    end

end

            