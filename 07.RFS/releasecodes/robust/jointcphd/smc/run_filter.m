function est = run_filter(model,meas)

% This is the MATLAB code for the lambda-pD-CPHD filter proposed in
% R. Mahler, B.-T. Vo and B.-N. Vo, "CPHD Filtering with unknown clutter rate and detection profile," IEEE Transactions on Signal Processing, Vol. 59, No. 8, pp. 3497-3513, 2011.
% http://ba-ngu.vo-au.com/vo/MVV_PHDrobust.pdf

% There are three versions of filters for this paper
% 1) Lambda-CPHD, 2) pD-CPHD, and 3) Lambda-pD-CPHD as sequentially described in the paper.
% This is the code for Lambda-pD-CPHD

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
est.L= zeros(meas.K,1);
est.pD= cell(meas.K,1);

%filter parameters
filter.J_max= 100000;                                         %total number of particles
filter.J_target= 1000;                                        %generated number of particles per expected target
filter.J_birth= model.L_birth*filter.J_target;                %generated number of particles from birth intensity
filter.Jc_birth= model.Lc_birth*filter.J_target;

filter.beta_factor = 1.1;           %beta factor


filter.N_max= 100;                   %maximum cardinality number (for cardinality distribution)

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%target initialization
w_update= eps;
m_init= [0.1;0;0.1;0;0.01];
P_init= diag([100 10 100 10 1]).^2;
u_update = 1;
v_update = 1;
x_update= [betarnd(u_update,v_update,1); m_init+chol(P_init)*randn(model.x_dim-1,1)];
l_update= 1;

%clutter initialization
P_D_init = 0.5;
Nc_init = round((size(meas.Z{1},2)-compute_pD(model,model.m_birth)*model.w_birth)/P_D_init);
cdn_update = [zeros(Nc_init-1,1);1;zeros(filter.N_max-Nc_init+1,1)];
Nc_update= Nc_init;


%recursive filtering
for k=1:meas.K
    %---intensity prediction 
    idx1= find(l_update==1);
    idx0= find(l_update==0);
    
    pS_vals= compute_pS_tg(model,x_update(:,idx1)); pS_vals= pS_vals(:);
    pS_vals_clt= compute_pS_clt(model,x_update(:,idx0)); pS_vals_clt= pS_vals_clt(:);

    x_predict= cat(2,gen_newstate_tg(model,x_update(:,idx1)), gen_newstate_clt(model,x_update(:,idx0)));
    w_temp = cat(1,pS_vals, pS_vals_clt);
    w_predict= w_temp.*w_update;
    l_predict= l_update;                                                                                                 

    for t=1:model.L_birth
        x_birth_temp1= repmat(model.m_birth(:,t), [1, filter.J_target])+model.B_birth(:,:,t)*randn(model.x_dim-1,filter.J_target);
        x_predict= [x_predict [betarnd(model.u_b(t),model.v_b(t),1,filter.J_target); x_birth_temp1]];                                                   %append birth particles
    end
    w_predict= cat(1,w_predict,sum(model.w_birth)*ones(filter.J_birth,1)/filter.J_birth);                                                               %append birth weights
    l_predict= cat(1,l_predict,ones(filter.J_birth,1));
    
    for t=1:model.Lc_birth
        x_birth_temp0= repmat(model.mc_birth(:,t), [1, filter.J_target])+model.Bc_birth(:,:,t)*randn(model.x_dim-1,filter.J_target);
        x_predict= [x_predict [betarnd(model.uc_b(t),model.vc_b(t),1,filter.J_target); x_birth_temp0]];                                                   %append birth particles
    end
    w_predict= cat(1,w_predict,sum(model.wc_birth)*ones(filter.Jc_birth,1)/filter.Jc_birth);                                                              %append birth weights
    l_predict= cat(1,l_predict, zeros(filter.Jc_birth,1));
    
    
    %---cardinality prediction 
    %surviving cardinality distribution
    survive_cdn_predict = zeros(filter.N_max+1,1);
    if ~isempty(idx0)
        survival_factor= sum(w_update(idx1).*pS_vals)/(sum(w_update(idx1))+sum(w_update(idx0))) + sum(w_update(idx0).*pS_vals_clt)/(sum(w_update(idx1))+sum(w_update(idx0)));
    else
        survival_factor= sum(w_update.*pS_vals)/(sum(w_update)+Nc_update);
    end
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
    for n=0:filter.N_max
        idxn=n+1;
        terms= zeros(filter.N_max+1,1);
        for j=0:n
            idxj= j+1;
            terms(idxj)= exp(-(sum(model.w_birth)+sum(model.wc_birth))+(n-j)*log(sum(model.w_birth)+sum(model.wc_birth))-sum(log(1:n-j)))*survive_cdn_predict(idxj);
        end
        cdn_predict(idxn) = sum(terms);
    end
    

    %normalize predicted cardinality distribution
    cdn_predict = cdn_predict/sum(cdn_predict);
    
    %--reorder track variables according to tags
    idx1= find(l_predict==1); 
    idx0= find(l_predict==0); 
    x_predict= [x_predict(:,idx1) x_predict(:,idx0)];
    l_predict= [l_predict(idx1);l_predict(idx0)];
    w_predict= [w_predict(idx1); w_predict(idx0)];
    
    
    %---intensity update
    %number of measurements
    m= size(meas.Z{k},2);    
    
    idx1= find(l_predict==1);
    idx0= find(l_predict==0);
    pD_vals= x_predict(1,idx1); pD_vals= pD_vals(:); qD_vals= 1-pD_vals;
    pD_vals_c= x_predict(1,idx0); pD_vals_c= pD_vals_c(:); qD_vals_c= 1-pD_vals_c;    
    
    %pre calculation for likelihood values
    if m~=0
        meas_likelihood1= zeros(length(w_predict(idx1)),m);
        meas_likelihood0= zeros(length(w_predict(idx0)),m);
        for ell=1:m
            meas_likelihood1(:,ell)= compute_likelihood_tg(model,meas.Z{k}(:,ell),x_predict(:,idx1))';
            meas_likelihood0(:,ell)= compute_likelihood_clt(model,meas.Z{k}(:,ell),x_predict(:,idx0))';
        end
    end
    %There is no elementary symmetrci functions in lambda-type chpd filter
    missed_factor= sum(w_predict(idx1).*qD_vals)/(sum(w_predict(idx1))+sum(w_predict(idx0))) + sum(w_predict(idx0).*qD_vals_c)/(sum(w_predict(idx1))+sum(w_predict(idx0)));
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

    w_update1 = 1/(sum(w_predict(idx1))+sum(w_predict(idx0)))*(qD_vals)*(Upsilon1'*cdn_predict)/(Upsilon0'*cdn_predict).*w_predict(idx1);
    w_update0= 1/(sum(w_predict(idx1))+sum(w_predict(idx0)))*(qD_vals_c)*(Upsilon1'*cdn_predict)/(Upsilon0'*cdn_predict).*w_predict(idx0);
    
    x_update1 = x_predict(:,idx1);
    x_update0= x_predict(:,idx0);
    l_update1= ones(size(w_update1,1),1); 
    l_update0= zeros(size(w_update0,1),1);
    
    % detection terms
    
    clut_w_term= zeros(m,1); w_cumsum= sum(w_update1);
    for ell=1:m   
       K_t= sum(pD_vals_c.*w_predict(idx0).*meas_likelihood0(:,ell));
       w_t = pD_vals.*w_predict(idx1).*meas_likelihood1(:,ell);        
       C= K_t + sum(w_t);
       w_t= w_t/C; w_cumsum= w_cumsum + sum(w_t);
       clut_w_term(ell)= K_t/C;
       w_update1= [w_update1; w_t];
       x_update1= [x_update1 x_predict(:,idx1)]; 
       l_update1= [l_update1; ones(size(w_t,1),1)];
    end
    w_update1= w_cumsum*w_update1/sum(w_update1);
    
    wc_cumsum= sum(w_update0);
    for ell=1:m   
       K_t= sum(pD_vals.*w_predict(idx1).*meas_likelihood1(:,ell));
       wc_t = pD_vals_c.*w_predict(idx0).*meas_likelihood0(:,ell);        
       C= K_t + sum(wc_t);
       wc_t= wc_t/C; wc_cumsum= wc_cumsum + sum(wc_t);
       w_update0= [w_update0; wc_t];
       x_update0= [x_update0 x_predict(:,idx0)]; 
       l_update0= [l_update0; ones(size(wc_t,1),1)];      
    end
    w_update0= wc_cumsum*w_update0/sum(w_update0);
    Nc_update= sum(w_predict(idx0))/(sum(w_predict)+sum(w_predict(idx0)))*(qD_vals_c'*w_predict(idx0))*(Upsilon1'*cdn_predict)/(Upsilon0'*cdn_predict) + sum(clut_w_term);

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
    
    
    %---resampling
    J_rsp= min(ceil(sum(w_update1)*filter.J_target),filter.J_max);
    idx= randsample(length(w_update1),J_rsp,true,w_update1); %idx= resample(w_update1,J_rsp);
    w_update1= sum(w_update1)*ones(J_rsp,1)/J_rsp;
    x_update1= x_update1(:,idx);
    
    J_rsp= min(ceil(sum(w_update0)*filter.J_target),filter.J_max);
    idx= randsample(length(w_update0),J_rsp,true,w_update0); %idx= resample(w_update0,J_rsp);
    w_update0= sum(w_update0)*ones(J_rsp,1)/J_rsp;
    x_update0= x_update0(:,idx);
    
    x_update= [x_update1 x_update0];
    w_update= [w_update1; w_update0];      w_posterior= w_update;
    l_update= [ones(size(w_update1,1),1);zeros(size(w_update0,1),1)];
    
    est.L(k) = x_update0(1,:)*w_update0(:);

    pD_tmp = [];
    %--- state extraction   
    if sum(w_update1) > .5
        [x_c,I_c]= our_kmeans(x_update1,w_update1,1);
        est.N(k)= 0;
        for j=1:size(x_c,2);
            if sum(w_update(I_c{j})) > .5,
                pD_tmp = [pD_tmp; x_update(1,I_c{j})*w_update(I_c{j})];
                est.X{k}= [ est.X{k} x_c(:,j) ];
                est.N(k)= est.N(k)+1;
            end
        end
    else
        est.N(k)= 0; est.X{k}= [];
    end
    est.pD{k} = w_update1(:)'*x_update1(1,:)'/sum(w_update1); %max(pD_tmp);
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #avg pD=' num2str(est.pD{k},4),...
         ' #eap lambda=' num2str(est.L(k),4),...
         ' #eap target=' num2str(sum(w_update1),4),...
         ' #est card=' num2str(est.N(k),4),...
         ' Neff_updt= ',num2str(round(1/sum((w_posterior/sum(w_posterior)).^2)))...
         ' Neff_rsmp= ',num2str(round(1/sum((w_update/sum(w_update)).^2)))   ]);
    end

end

            