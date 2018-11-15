function est = run_filter(model,meas)

% This is the MATLAB code for the CBMeMBer filter proposed in
% (without track labelling)
% B.-T. Vo, B.-N. Vo, and A. Cantoni, "The Cardinality Balanced Multi-target Multi-Bernoulli filter and its implementations," IEEE Trans. Signal Processing, Vol. 57, No. 2, pp. 409–423, 2009. 
% http://ba-ngu.vo-au.com/vo/VVCmemberSP09.pdf
% ---BibTeX entry
% @ARTICLE{CBMEMBER,
% author={B.-T. Vo and B.-N. Vo and A. Cantoni},
% journal={IEEE Transactions on Signal Processing},
% title={The Cardinality Balanced Multi-Target Multi-Bernoulli Filter and Its Implementations},
% year={2009},
% month={Feb},
% volume={57},
% number={2},
% pages={409-423}} 
%---

%=== Setup

%output variables
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);
est.L= cell(meas.K,1);

%filter parameters
filter.T_max= 100;                  %maximum number of tracks
filter.track_threshold= 1e-3;       %threshold to prune tracks

filter.L_max= 100;                  %limit on number of Gaussians in each track
filter.elim_threshold= 1e-5;        %pruning threshold for Gaussians in each track
filter.merge_threshold= 4;          %merging threshold for Gaussians in each track

filter.P_G= 0.9999999;                           %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 0;                             %gating on or off 1/0

filter.ukf_alpha= 1;                %scale parameter for UKF - choose alpha=1 ensuring lambda=beta and offset of first cov weight is beta for numerical stability
filter.ukf_beta= 2;                 %scale parameter for UKF
filter.ukf_kappa= 2;                %scale parameter for UKF (alpha=1 preferred for stability, giving lambda=1, offset of beta for first cov weight)
filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
T_update= 1;
r_update(1)= 0.001;
w_update{1}(1)= 1;
m_update{1}(:,1)= [0.1;0;0.1;0;0.01];
P_update{1}(:,:,1)= diag([100 10 100 10 1]).^2;
L_update(1)= 1;

%recursive filtering
for k=1:meas.K
    %---prediction 
    T_predict= T_update+model.T_birth;                                                                              %total number of tracks/components
    L_predict= zeros(T_predict,1);                                                                                  %total number of Gaussians in each track/component
    r_predict= zeros(T_predict,1);                                                                                  %existence probability for tracks/components
    m_predict= cell(T_predict,1);                                                                                   %means of Gaussians in each track/component
    P_predict= cell(T_predict,1);                                                                                   %covs of Gaussians in each track/component
    w_predict= cell(T_predict,1);                                                                                   %weights of Gaussians in each track/component
    
    %prediction for surviving tracks
    for t=1:T_update
        L_predict(t)= L_update(t);                                                                                                                              %surviving number of Gaussians in current track
        r_predict(t)= model.P_S*r_update(t);                                                                                                                    %surviving existence probability
        [m_predict{t},P_predict{t}] = ukf_predict_multiple(model,m_update{t},P_update{t},filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);                    %surviving components
        w_predict{t}= w_update{t};                                                                                                                              %surviving weights
    end                                                  

    %insertion of birth tracks
    offset= T_update;
    for t=1:model.T_birth
        L_predict(offset+t)= model.L_birth(1);                                                                                                                  %append birth Gaussian count
        r_predict(offset+t)= model.r_birth(t);                                                                                                                  %append birth probabilities
        m_predict{offset+t}= model.m_birth{t}; P_predict{offset+t}= model.P_birth{t};                                                                           %append birth components
        w_predict{offset+t}= model.w_birth{t};                                                                                                                  %append birth weights 
    end                     

    r_predict= limit_range(r_predict);                                                                              %limit range of 0<r<1 for numerical stability

    %---construction of pseudo PHD for update
    L_pseudo= sum(L_predict);                                           %number of Gaussians in pseudo-PHD
    m_pseudo= zeros(model.x_dim,L_pseudo);                              %means of Gaussians in pseudo-PHD
    P_pseudo= zeros(model.x_dim,model.x_dim,L_pseudo);                  %covs of Gaussians in pseudo-PHD
    w_pseudo= zeros(L_pseudo,1);                                        %weights of Gaussians in pseudo-PHD
    w_pseudo1= zeros(L_pseudo,1);                                       %alt weight (1) of Gaussians in pseudo-PHD - used in CB-MeMBer update later
    w_pseudo2= zeros(L_pseudo,1);                                       %alt weight (2) of Gaussians in pseudo-PHD - used in CB-MeMBer update later
    
    start_pt= 1;
    for t=1:T_predict                  
        end_pt= start_pt+L_predict(t)-1;
        m_pseudo(:,start_pt:end_pt)= m_predict{t};
        P_pseudo(:,:,start_pt:end_pt)= P_predict{t};
        w_pseudo(start_pt:end_pt) = r_predict(t)/(1-r_predict(t))*w_predict{t};
        w_pseudo1(start_pt:end_pt)= r_predict(t)/(1-r_predict(t)*model.P_D)*w_predict{t};
        w_pseudo2(start_pt:end_pt)= r_predict(t)*(1-r_predict(t))/((1-r_predict(t)*model.P_D)^2)*w_predict{t};
        start_pt= end_pt+1;
    end
    
    %---gating
    if filter.gate_flag
        meas.Z{k}= gate_meas_ukf(meas.Z{k},filter.gamma,model,m_pseudo,P_pseudo,filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);       
    end
        
    %---update
    %number of measurements
    m= size(meas.Z{k},2);
    
    T_update= T_predict+m;                                                                                        %total number of tracks/components
    L_update= zeros(T_update,1);                                                                                  %total number of Gaussians in each track/component
    r_update= zeros(T_update,1);                                                                                  %existence probability for tracks/components
    m_update= cell(T_update,1);                                                                                   %means of Gaussians in each track/component
    P_update= cell(T_update,1);                                                                                   %covs of Gaussians in each track/component
    w_update= cell(T_update,1);                                                                                   %weights of Gaussians in each track/component
    
    %legacy tracks
    L_update= L_predict;
    r_update= r_predict.*((1-model.P_D)./(1-r_predict*model.P_D)); 
    m_update= m_predict;
    P_update= P_predict;
    w_update= w_predict;
    
    %measurement updated tracks
    if m~=0
        offset= T_predict;
        [qz_temp,m_temp,P_temp] = ukf_update_multiple(meas.Z{k},model,m_pseudo,P_pseudo,filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);
        for ell=1:m                                        
            r_temp= (model.P_D*w_pseudo2(:)'*qz_temp(:,ell))/(model.lambda_c*model.pdf_c+model.P_D*w_pseudo1(:)'*qz_temp(:,ell));
            w_temp= model.P_D*w_pseudo(:).*qz_temp(:,ell); w_temp= w_temp/sum(w_temp);
            L_temp= length(w_temp);
            
            L_update(offset+ell)= L_temp; 
            r_update(offset+ell)= r_temp;
            m_update{offset+ell}= m_temp(:,:,ell);
            P_update{offset+ell}= P_temp;
            w_update{offset+ell}= w_temp;    
        end
        
    end
    
    T_update= T_predict+m;
    r_update= limit_range(r_update);    
            
    %---track/component management
    T_posterior= T_update;
    
    %tracks: pruning, capping
    [r_update,w_update,m_update,P_update,L_update]= prune_tracks(r_update,w_update,m_update,P_update,L_update,filter.track_threshold);  T_prune= length(r_update);
    [r_update,w_update,m_update,P_update,L_update]= cap_tracks(r_update,w_update,m_update,P_update,L_update,filter.T_max);              T_cap  = length(r_update);  
    
    T_update= T_cap;
    
    %(within tracks) pruning,merging,capping
    for t=1:T_update
        [w_update{t},m_update{t},P_update{t}]= gaus_prune(w_update{t},m_update{t},P_update{t},filter.elim_threshold);    L_prune(t)= length(w_update{t});
        [w_update{t},m_update{t},P_update{t}]= gaus_merge(w_update{t},m_update{t},P_update{t},filter.merge_threshold);   L_merge(t)= length(w_update{t});
        [w_update{t},m_update{t},P_update{t}]= gaus_cap(w_update{t},m_update{t},P_update{t},filter.L_max);               L_cap(t)  = length(w_update{t});
    end
    
    L_update= L_cap;
    
    %--- state extraction
    cdn_update= prod(1-r_update)*esf(r_update./(1-r_update));
    [~,idx_max_cdn] = max(cdn_update);
    map_cdn = idx_max_cdn-1;
    est.N(k) = min(length(r_update),map_cdn);
    [~,idx_max_tgt]= sort(-r_update);
    for t=1:est.N(k)
         [~,idx]= max(w_update{idx_max_tgt(t)});
        est.X{k} = [est.X{k} m_update{idx_max_tgt(t)}(:,idx)];
        est.L{k}= [];
    end
  
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #eap cdn=' num2str(sum(r_update)),...
         ' #var cdn=' num2str(sum(r_update.*(1-r_update)),4),...
         ' #est card=' num2str(est.N(k),4),...
         ' #trax updt=' num2str(T_posterior,4),...
         ' #trax elim=' num2str(T_prune,4),...
         ' #trax filt=' num2str(T_cap,4),...
         ' #gaus filt=',num2str(sum(L_cap))   ]);
    end

end

function clipped_r= limit_range(r)

r(r>0.999)=0.999;
r(r<0.001)=0.001;
clipped_r= r;

function [new_r,new_w,new_m,new_P,new_L]= prune_tracks(r,w,m,P,L,track_threshold)

idx= find(r>track_threshold);

new_r= zeros(length(idx));
new_m= cell(length(idx),1);
new_P= cell(length(idx),1);
new_w= cell(length(idx),1);
new_L= zeros(length(idx),1);

new_r= r(idx);
for i=1:length(idx)
    new_m{i}= m{idx(i)};
    new_P{i}= P{idx(i)};
    new_w{i}= w{idx(i)};
end
new_L= L(idx);

function [new_r,new_w,new_m,new_P,new_L]= cap_tracks(r,w,m,P,L,T_max)

if length(r) > T_max
    
    [~,idx]= sort(-r);
    
    new_r= zeros(length(idx));
    new_m= cell(length(idx),1);
    new_P= cell(length(idx),1);
    new_w= cell(length(idx),1);
    new_L= zeros(length(idx),1);
    
    new_r= r(idx(1:T_max));
    for i=1:T_max
        new_m{i}= m{idx(i)};
        new_P{i}= P{idx(i)};
        new_w{i}= w{idx(i)};
    end
    new_L= L(idx(1:T_max));
    
else
    new_r= r;
    new_m= m;
    new_P= P;
    new_w= w;
    new_L= L;
end




            