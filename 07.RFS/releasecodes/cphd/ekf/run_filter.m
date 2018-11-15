function est = run_filter(model,meas)

% This is the MATLAB code for the GMCPHD filter proposed in
% (assuming Poisson clutter)
% B.-T. Vo, B.-N. Vo and A. Cantoni, "Analytic implementations of the Cardinalized Probability Hypothesis Density Filter," IEEE Trans Signal Processing, Vol. 55,  No. 7, part 2,  pp. 3553-3567, 2007.
% http://ba-ngu.vo-au.com/vo/VVC_CPHD_SP07.pdf
% ---BibTeX entry
% @ARTICLE{GMCPHD,
% author={B.-T. Vo and B.-N. Vo and A. Cantoni},
% journal={IEEE Transactions on Signal Processing},
% title={Analytic Implementations of the Cardinalized Probability Hypothesis Density Filter},
% year={2007},
% month={July},
% volume={55},
% number={7},
% pages={3553-3567}} 
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

filter.N_max= 20;                   %maximum cardinality number (for cardinality distribution)

filter.P_G= 0.9999999;                           %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 1;                             %gating on or off 1/0

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
w_update(1)= eps;
m_update(:,1)= [0.1;0;0.1;0;0];
P_update(:,:,1)= diag([1 1 1 1 0.1]).^2;
L_update = 1;
cdn_update= [1; zeros(filter.N_max,1)];

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
    
    %---gating
    if filter.gate_flag
        meas.Z{k}= gate_meas_ekf(meas.Z{k},filter.gamma,model,m_predict,P_predict);        
    end
        
    %---intensity update
    %number of measurements
    m= size(meas.Z{k},2);
    
    %pre calculation for Kalman update parameters
    if m~=0
        [qz_temp,m_temp,P_temp] = ekf_update_multiple(meas.Z{k},model,m_predict,P_predict);
    end
    
    %pre calculation for elementary symmetric functions
    XI_vals = zeros(m,1);                        %arguments to esf
    for ell=1:m
       XI_vals(ell) = model.P_D*w_predict'*qz_temp(:,ell)/model.pdf_c;
    end
    
    esfvals_E = esf(XI_vals);                   %calculate esf for entire observation set
    esfvals_D = zeros(m,m);                     %calculate esf with each observation index removed one-by-one
    for ell=1:m
        esfvals_D(:,ell) = esf([XI_vals(1:ell-1);XI_vals(ell+1:m)]);
    end
    
    %pre calculation for upsilons
    upsilon0_E = zeros(filter.N_max+1,1);
    upsilon1_E = zeros(filter.N_max+1,1);
    upsilon1_D = zeros(filter.N_max+1,m);
    
    for n=0:filter.N_max
        idxn= n+1;
        
        terms0_E= zeros(min(m,n)+1,1);  %calculate upsilon0_E(idxn)
        for j=0:min(m,n)
            idxj= j+1;
            terms0_E(idxj) = exp(-model.lambda_c+(m-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-j))+(n-j)*log(model.Q_D)-j*log(sum(w_predict)))*esfvals_E(idxj);
        end
        upsilon0_E(idxn)= sum(terms0_E);
        
        terms1_E= zeros(min(m,n)+1,1);  %calculate upsilon1_E(idxn)
        for j=0:min(m,n)
            idxj= j+1;
            if n>=j+1
                terms1_E(idxj) = exp(-model.lambda_c+(m-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-(j+1)))+(n-(j+1))*log(model.Q_D)-(j+1)*log(sum(w_predict)))*esfvals_E(idxj);
            end
        end
        upsilon1_E(idxn)= sum(terms1_E);
        
        if m~= 0                        %calculate upsilon1_D(idxn,:) if m>0
            terms1_D= zeros(min((m-1),n)+1,m);
            for ell=1:m
                for j=0:min((m-1),n)
                    idxj= j+1;
                    if n>=j+1
                        terms1_D(idxj,ell) = exp(-model.lambda_c+((m-1)-j)*log(model.lambda_c)+sum(log(1:n))-sum(log(1:n-(j+1)))+(n-(j+1))*log(model.Q_D)-(j+1)*log(sum(w_predict)))*esfvals_D(idxj,ell);
                    end
                end
            end
            upsilon1_D(idxn,:)= sum(terms1_D,1);
        end
    end
    

    %missed detection term 
    w_update = (upsilon1_E'*cdn_predict)/(upsilon0_E'*cdn_predict)*model.Q_D*w_predict;
    m_update = m_predict;
    P_update = P_predict;
    
    if m~=0
        %m detection terms 
        %Kalman update precalculated and stored
        for ell=1:m
            w_temp = (upsilon1_D(:,ell)'*cdn_predict)/(upsilon0_E'*cdn_predict)*model.P_D.*qz_temp(:,ell)/model.pdf_c.*w_predict;
            w_update = cat(1,w_update,w_temp);
            m_update = cat(2,m_update,m_temp(:,:,ell));
            P_update = cat(3,P_update,P_temp);
        end
    end    
    
    %---cardinality update
    cdn_update= upsilon0_E.*cdn_predict;
    cdn_update= cdn_update/sum(cdn_update);
            
    %---mixture management
    L_posterior= length(w_update);
    
    %pruning, merging, capping
    [w_update,m_update,P_update]= gaus_prune(w_update,m_update,P_update,filter.elim_threshold);    L_prune= length(w_update);
    [w_update,m_update,P_update]= gaus_merge(w_update,m_update,P_update,filter.merge_threshold);   L_merge= length(w_update);
    [w_update,m_update,P_update]= gaus_cap(w_update,m_update,P_update,filter.L_max);               L_cap  = length(w_update);
    
    L_update= L_cap;
    
    %--- state extraction
    [~,idx_max_cdn] = max(cdn_update);
    map_cdn = idx_max_cdn-1;
    est.N(k) = min(length(w_update),map_cdn);
    [~,idx_m_srt]= sort(-w_update);
    est.X{k} = m_update(:,idx_m_srt(1:est.N(k)));
    
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #tot phd=' num2str(sum(w_update),4),...
         ' #eap cdn=' num2str([0:filter.N_max]*cdn_update,4),...
         ' #var cdn=' num2str([0:filter.N_max].^2*cdn_update-([0:filter.N_max]*cdn_update)^2,4),...
         ' #est card=' num2str(est.N(k),4),...
         ' #gaus orig=',num2str(L_posterior),...
         ' #gaus elim=',num2str(L_prune), ...
         ' #gaus merg=',num2str(L_merge)   ]);
    end

end

            