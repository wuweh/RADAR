function est = run_filter(model,meas)

% This is the MATLAB code for the lambda-CPHD filter proposed in
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

%=== Setup

%output variables
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);
est.L= cell(meas.K,1);

%filter parameters
filter.L_max= 100;                  %limit on number of Gaussians
filter.elim_threshold= 1e-5;        %pruning threshold
filter.merge_threshold= 4;          %merging threshold

filter.N_max= 300;                   %maximum cardinality number (for cardinality distribution)

filter.P_G= 0.9999999;                           %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 1;                             %gating on or off 1/0

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initialisation of targets
w_update(1)= eps;
m_update(:,1)= [0.1;0;0.1;0;0];
P_update(:,:,1)= diag([1 1 1 1 0.1]).^2;
L_update = 1;

%initialisation of clutter
Nc_init = round((size(meas.Z{1},2)-model.P_D*sum(model.w_birth))/model.clutter_P_D);
cdn_update = [zeros(Nc_init-1,1);1;zeros(filter.N_max-Nc_init+1,1)];
Nc_update = Nc_init;

%recursive filtering
for k=1:meas.K
    %---intensity prediction 
    [m_predict,P_predict] = ekf_predict_multiple(model,m_update,P_update);                       %surviving components
    w_predict= model.P_S*w_update;                                                                  %surviving weights

    m_predict= cat(2,model.m_birth,m_predict); P_predict=cat(3,model.P_birth,P_predict);            %append birth components
    w_predict= cat(1,model.w_birth,w_predict);                                                      %append birth weights
                                                 
    L_predict= model.L_birth+L_update;                                                              %number of predicted components
    
    %---cardinality prediction 
    %surviving cardinality distribution
    survive_cdn_predict = zeros(filter.N_max+1,1);
    survival_factor= sum(w_update)/(sum(w_update)+Nc_update)*model.P_S + Nc_update/(sum(w_update)+Nc_update)*model.clutter_P_S;
    if isnan(survival_factor), survival_factor=0; end %catch the degernerate zero case
    for j=0:filter.N_max
        idxj=j+1;
        terms= zeros(filter.N_max+1,1);
        for ell=j:filter.N_max
            idxl= ell+1;
            terms(idxl) = exp(sum(log(1:ell))-sum(log(1:j))-sum(log(1:ell-j))+j*log(survival_factor)+(ell-j)*log(1-survival_factor))*cdn_update(idxl);
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
    %normalize predicted cardinality distribution
    cdn_predict = cdn_predict/sum(cdn_predict);
    Nc_predict = model.lambda_cb + model.clutter_P_S * Nc_update;
    
           
    %---intensity update
    %number of measurements
    m= size(meas.Z{k},2);
    
    %pre calculation for Kalman update parameters
    if m~=0
        [qz_temp,m_temp,P_temp] = ekf_update_multiple(meas.Z{k},model,m_predict,P_predict);
    end
    
    %There is no elementary symmetric functions in lambda-CPHD filter
    %pre calculation for Upsilon0 and Upsilon1
    
    missed_factor= sum(w_predict)/(sum(w_predict)+Nc_predict)*(1-model.P_D) + Nc_predict/(sum(w_predict)+Nc_predict)*(1-model.clutter_P_D);
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
 
    % misdetection term
    m_update = m_predict;
    P_update = P_predict;
    w_update = 1/(sum(w_predict)+Nc_predict)*(1-model.P_D)*(Upsilon1'*cdn_predict)/(Upsilon0'*cdn_predict)*w_predict;
     
    % detection terms
    clut_w_term= zeros(m,1); w_cumsum= sum(w_update);
    start_pt= L_predict+1;
    for ell=1:m
        
       %calculate weights first 
       K_t= Nc_predict*model.clutter_P_D/prod(model.range_c(:,2)-model.range_c(:,1));
       w_t = model.P_D*w_predict.*qz_temp(:,ell);
       C= K_t + sum(w_t);
       w_t= w_t/C; w_cumsum= w_cumsum + sum(w_t);
       clut_w_term(ell)= K_t/C;
        
       %find terms with weights that won't be truncated at the end
       idx= find( w_t > filter.elim_threshold ); %idx= find( w_t > 0 ); %if you want to do them all
       end_pt= start_pt-1 + length(idx);
       w_update(start_pt:end_pt)= w_t(idx);
        
       %update for these terms
       for j=1:length(idx)
           m_update(:,start_pt-1+j) = m_temp(:,idx(j),ell);
       end
       P_update(:,:,start_pt:end_pt)= P_temp(:,:,idx);
       start_pt= end_pt+ 1;
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
    
    % mixture management
    L_posterior = length(w_update);
    
    %pruning, merging, capping
    [w_update,m_update,P_update]= gaus_prune(w_update,m_update,P_update,filter.elim_threshold);    L_prune= length(w_update);
    [w_update,m_update,P_update]= gaus_merge(w_update,m_update,P_update,filter.merge_threshold);   L_merge= length(w_update);
    [w_update,m_update,P_update]= gaus_cap(w_update,m_update,P_update,filter.L_max);               L_cap  = length(w_update);
    
    L_update= L_cap;
   
    
    %--- state extraction
    est.N_soft(k) = sum(w_update);
    est.N(k) = min(length(w_update),round(est.N_soft(k)));
    [~,idx_m_srt]= sort(-w_update);
    est.X{k} = m_update(:,idx_m_srt(1:est.N(k)));
    est.L{k} = Nc_update*model.clutter_P_D;
    
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #eap lambda=' num2str(est.L{k},4),...
         ' #eap target=' num2str(sum(w_update),4),...
         ' #est card=' num2str(est.N(k),4),...
         ' #gaus orig=',num2str(L_posterior),...
         ' #gaus elim=',num2str(L_prune), ...
         ' #gaus merg=',num2str(L_merge)   ]);
    end

end

            