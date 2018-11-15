function est = run_filter(model,meas)

% This is the MATLAB code for the pD-CPHD filter proposed in
% R. Mahler, B.-T. Vo and B.-N. Vo, "CPHD Filtering with unknown clutter rate and detection profile," IEEE Transactions on Signal Processing, Vol. 59, No. 8, pp. 3497-3513, 2011.
% http://ba-ngu.vo-au.com/vo/MVV_PHDrobust.pdf

% There are three versions of filters for this paper
% 1) Lambda-CPHD, 2) pD-CPHD, and 3) Lambda-pD-CPHD as sequentially described in the paper.
% This is the code for pD-CPHD

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
est.pD= cell(meas.K,1);

%filter parameters
filter.L_max= 1000;                 %limit on number of Gaussians
filter.elim_threshold= 1e-5;        %pruning threshold
filter.merge_threshold= 4;          %merging threshold
filter.beta_factor = 1.1;           %beta factor

filter.N_max= 20;                   %maximum cardinality number (for cardinality distribution)

filter.P_G= 0.999;                               %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 1;                             %gating on or off 1/0

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
w_update(1)= eps;
m_update(:,1)= [0.1;0;0.1;0];
P_update(:,:,1)= diag([1 1 1 1]).^2;
u_update = 1;
v_update = 1;
L_update = 1;
cdn_update= [1; zeros(filter.N_max,1)];

%recursive filtering
for k=1:meas.K
    %---intensity prediction 
    [m_predict,P_predict] = kalman_predict_multiple(model,m_update,P_update);                       %surviving components
    w_predict= model.P_S*w_update;                                                                  %surviving weights
    
    %--- prediction of beta parameters
    u_predict = zeros(L_update,1); v_predict = zeros(L_update,1);
    for j=1:L_update
        beta_avg_tmp= u_update(j)/(u_update(j)+v_update(j)); 
        beta_var_tmp= (u_update(j)*v_update(j))/((u_update(j)+v_update(j))^2*(u_update(j)+v_update(j)+1)); 
        beta_tht_tmp= (beta_avg_tmp*(1-beta_avg_tmp))/min(filter.beta_factor*beta_var_tmp,beta_avg_tmp*(1-beta_avg_tmp))- 1;
        u_predict(j)= beta_tht_tmp*beta_avg_tmp;  v_predict(j)= beta_tht_tmp*(1-beta_avg_tmp);
    end   
    
    m_predict= cat(2,model.m_birth,m_predict); P_predict=cat(3,model.P_birth,P_predict);            %append birth components
    w_predict= cat(1,model.w_birth,w_predict);                                                      %append birth weights
    u_predict = cat(1,model.u_b,u_predict);
    v_predict = cat(1,model.v_b,v_predict);
    
    L_predict= model.L_birth+L_update;                                                              %number of predicted components
    
    %---cardinality prediction 
    %surviving cardinality distribution
    survive_cdn_predict = zeros(filter.N_max+1,1);
    for j=0:filter.N_max
        idxj=j+1;
        terms= zeros(filter.N_max+1,1);
        for ell=j:filter.N_max
            idxl= ell+1;
            terms(idxl) = exp(sum(log(1:ell))-sum(log(1:j))-sum(log(1:ell-j))+j*log(model.P_S)+(ell-j)*log(model.Q_S))*cdn_update(idxl);
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
            terms(idxj)= exp(-sum(model.w_birth)+(n-j)*log(sum(model.w_birth))-sum(log(1:n-j)))*survive_cdn_predict(idxj);
        end
        cdn_predict(idxn) = sum(terms);
    end

    %normalize predicted cardinality distribution
    cdn_predict = cdn_predict/sum(cdn_predict);
    
        
    %---intensity update
    %number of measurements
    m= size(meas.Z{k},2);
    
    %pre calculation for Kalman update parameters
    if m~=0
        [qz_temp,m_temp,P_temp] = kalman_update_multiple(meas.Z{k},model,m_predict,P_predict);
    end
    
    %pre calculation for elementary symmetric functions
    % mean Beta-PD values
    avg_BPD= u_predict./(u_predict+v_predict);
    % integral between beta-prior-likelihood
    intgv = zeros(m,1);
    for ell=1:m
        intgv(ell) = (w_predict'.*avg_BPD')*qz_temp(:,ell);
    end
    zvals = intgv./model.pdf_c;
    esfvals_E = esf(zvals); % calcuate esf for entire observation set
    esfvals_D = zeros(m,m);
    for ell=1:m
        esfvals_D(:,ell) = esf([zvals(1:ell-1);zvals(ell+1:m)]);
    end 
    
    upsilon0_E = zeros(filter.N_max+1,1);
    upsilon1_E = zeros(filter.N_max+1,1);
    upsilon1_D = zeros(filter.N_max+1,m);
    
    d_factor = sum(w_predict.*avg_BPD/sum(w_predict));
    for n=0:filter.N_max
        idxn= n+1;
        
        terms0_E= zeros(min(m,n)+1,1);  %calculate upsilon0_E(idxn)
        for j=0:min(m,n)
            idxj= j+1;
            terms0_E(idxj) = exp(-model.lambda_c+(m-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-j))+(n-j)*log(1-d_factor)-j*log(sum(w_predict)))*esfvals_E(idxj);
        end
        upsilon0_E(idxn)= sum(terms0_E);
        
        terms1_E= zeros(min(m,n)+1,1);  %calculate upsilon1_E(idxn)
        for j=0:min(m,n)
            idxj= j+1;
            if n>=j+1
                terms1_E(idxj) = exp(-model.lambda_c+(m-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-(j+1)))+(n-(j+1))*log(1-d_factor)-(j+1)*log(sum(w_predict)))*esfvals_E(idxj);
            end
        end
        upsilon1_E(idxn)= sum(terms1_E);
        
        if m~= 0                        %calculate upsilon1_D(idxn,:) if m>0
            terms1_D= zeros(min((m-1),n)+1,m);
            for ell=1:m
                for j=0:min((m-1),n)
                    idxj= j+1;
                    if n>=j+1
                        terms1_D(idxj,ell) = exp(-model.lambda_c+((m-1)-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-(j+1)))+(n-(j+1))*log(1-d_factor)-(j+1)*log(sum(w_predict)))*esfvals_D(idxj,ell);
                    end
                end
            end
            upsilon1_D(idxn,:)= sum(terms1_D,1);
        end
    end
    

    %missed detection term 
    w_update = (upsilon1_E'*cdn_predict)/(upsilon0_E'*cdn_predict)*(beta(u_predict,v_predict+1)./beta(u_predict,v_predict).*w_predict);
    m_update = m_predict;
    P_update = P_predict;
    u_update = u_predict; v_update = v_predict+1;
    
    if m~=0
        %m detection terms 
        %Kalman update precalculated and stored
        for ell=1:m
            u_temp = u_predict+1; v_temp = v_predict;
            w_temp = (upsilon1_D(:,ell)'*cdn_predict)/(upsilon0_E'*cdn_predict)*(beta(u_predict+1,v_predict)./beta(u_predict,v_predict)).*w_predict.*qz_temp(:,ell)./model.pdf_c;
            w_update = cat(1,w_update,w_temp);
            m_update = cat(2,m_update,m_temp(:,:,ell));
            P_update = cat(3,P_update,P_temp);
            u_update = cat(1,u_update, u_temp);
            v_update = cat(1,v_update, v_temp);
        end
    end    
    
    %---cardinality update
    cdn_update= upsilon0_E.*cdn_predict;
    cdn_update= cdn_update/sum(cdn_update);
            
    %---mixture management
    L_posterior= length(w_update);
    
    %pruning, merging, capping
    [w_update,m_update,P_update,u_update,v_update]= gaus_prune_1(w_update,m_update,P_update,u_update,v_update, filter.elim_threshold);    L_prune= length(w_update);
    L_merge= length(w_update); %merging disabled due to higher computational cost
    [w_update,m_update,P_update,u_update,v_update]= gaus_cap_1(w_update,m_update,P_update,u_update,v_update,filter.L_max);               L_cap  = length(w_update);
    
    L_update= L_cap;
    
    %--- state extraction
    [~,idx_max_cdn] = max(cdn_update);
    map_cdn = idx_max_cdn-1;
    est.N(k) = min(length(w_update),map_cdn);
    [~,idx_m_srt]= sort(-w_update);
    est.X{k} = m_update(:,idx_m_srt(1:est.N(k)));
    pD_tmp = (u_update(:)'./(u_update(:)'+v_update(:)'));
    est.pD{k}= w_update(:)'*pD_tmp'/sum(w_update);

    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #avg pD =' num2str(est.pD{k},4),...
         ' #eap target=' num2str(sum(w_update),4),...
         ' #est card=' num2str(est.N(k),4),...
         ' #gaus orig=',num2str(L_posterior),...
         ' #gaus elim=',num2str(L_prune), ...
         ' #gaus merg=',num2str(L_merge)   ]);

    end

end

            