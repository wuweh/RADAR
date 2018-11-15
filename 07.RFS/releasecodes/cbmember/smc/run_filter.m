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

filter.J_max= 1000;                  %max number of particles in each track (otherwise allocated proportion to existence probability)
filter.J_min= 300;                   %min number of particles in each track (otherwise allocated proportion to existence probability)

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
T_update= 1;
r_update(1)= 0.001;
w_update{1}= ones(filter.J_min,1)/filter.J_min;
m_init= [0.1;0;0.1;0;0.01]; P_init= diag([100 10 100 10 1]).^2;
x_update{1}= gen_gms(1,m_init,P_init,filter.J_min);
J_update(1)= filter.J_min;

%recursive filtering
for k=1:meas.K
    %---prediction 
    T_predict= T_update+model.T_birth;                                                                              %total number of tracks/components
    J_predict= zeros(T_predict,1);                                                                                  %total number of particles in each track/component
    r_predict= zeros(T_predict,1);                                                                                  %existence probability for tracks/components
    x_predict= cell(T_predict,1);                                                                                   %states of particles in each track/component
    w_predict= cell(T_predict,1);                                                                                   %weights of particles in each track/component
    
    %prediction for surviving tracks
    for t=1:T_update
        w_temp= compute_pS(model,x_update{t}).*w_update{t};
        J_predict(t)= J_update(t);                                                                                                                              %surviving number of particles in current track
        r_predict(t)= r_update(t)*sum(w_temp);                                                                                                                  %surviving existence probability
        x_predict{t}= gen_newstate_fn(model,x_update{t},'noise');                                                                                               %surviving particles
        w_predict{t}= w_temp/sum(w_temp);                                                                                                                       %surviving weights

    end                                                  

    %insertion of birth tracks
    offset= T_update;
    for t=1:model.T_birth
        J_predict(offset+t)= max(round(model.r_birth(t)*filter.J_max),filter.J_min);                                                                            %append birth particles count
        r_predict(offset+t)= model.r_birth(t);                                                                                                                  %append birth probabilities
        x_predict{offset+t}= gen_gms(model.w_birth{t},model.m_birth{t},model.P_birth{t},J_predict(offset+t));                                                   %append birth particles
        w_predict{offset+t}= ones(J_predict(offset+t),1)/J_predict(offset+t);                                                                                   %append birth weights
    end                     

    r_predict= limit_range(r_predict);                                                                              %limit range of 0<r<1 for numerical stability

    %---construction of pseudo PHD for update
    J_pseudo= sum(J_predict);                                           %number of particles in pseudo-PHD
    x_pseudo= zeros(model.x_dim,J_pseudo);                              %states of particles in pseudo-PHD
    w_pseudo= zeros(J_pseudo,1);                                        %weights of particles in pseudo-PHD
    w_pseudo1= zeros(J_pseudo,1);                                       %alt weight (1) of particles in pseudo-PHD - used in CB-MeMBer update later
    w_pseudo2= zeros(J_pseudo,1);                                       %alt weight (2) of particles in pseudo-PHD - used in CB-MeMBer update later
    
    start_pt= 1;
    for t=1:T_predict                  
        end_pt= start_pt+J_predict(t)-1;
        w_temp= compute_pD(model,x_predict{t}).*w_predict{t};
        x_pseudo(:,start_pt:end_pt)= x_predict{t};
        w_pseudo(start_pt:end_pt) = r_predict(t)/(1-r_predict(t))*w_predict{t};
        w_pseudo1(start_pt:end_pt)= r_predict(t)/(1-r_predict(t)*sum(w_temp))*w_predict{t};
        w_pseudo2(start_pt:end_pt)= r_predict(t)*(1-r_predict(t))/((1-r_predict(t)*sum(w_temp))^2)*w_predict{t};
        start_pt= end_pt+1;
    end
    
    %---update
    %number of measurements
    m= size(meas.Z{k},2);
    
    T_update= T_predict+m;                                                                                        %total number of tracks/components
    J_update= zeros(T_update,1);                                                                                  %total number of particles in each track/component
    r_update= zeros(T_update,1);                                                                                  %existence probability for tracks/components
    x_update= cell(T_update,1);                                                                                   %states of particles in each track/component
    w_update= cell(T_update,1);                                                                                   %weights of particles in each track/component
 
    %legacy tracks
    for t=1:T_predict  
        w_temp= w_predict{t}.*compute_pD(model,x_predict{t});
        J_update(t)= J_predict(t);
        r_update(t)= r_predict(t)*((1-sum(w_temp))/(1-r_predict(t)*sum(w_temp)));
        x_update{t}= x_predict{t};
        w_update{t}= w_predict{t}.*(1-compute_pD(model,x_predict{t}))/(1-sum(w_temp));
    end
    
    %measurement updated tracks
    if m~=0
        offset= T_predict;
        pD_vals= compute_pD(model,x_pseudo);
        for ell=1:m
            meas_likelihood= compute_likelihood(model,meas.Z{k}(:,ell),x_pseudo)';
            r_temp= sum(pD_vals.*w_pseudo2(:).*meas_likelihood)/(model.lambda_c*model.pdf_c+sum(pD_vals.*w_pseudo1(:).*meas_likelihood));
            w_temp= pD_vals.*w_pseudo.*meas_likelihood; w_temp= w_temp/sum(w_temp);
            J_temp= length(w_temp);
            
            J_update(offset+ell)= J_temp; 
            r_update(offset+ell)= r_temp;
            x_update{offset+ell}= x_pseudo;
            w_update{offset+ell}= w_temp; 
        end
        
    end

    r_update= limit_range(r_update);    
            
    %---track/component management
    T_posterior= T_update;
    
    %tracks: pruning, capping
    [r_update,w_update,x_update,J_update]= prune_tracks(r_update,w_update,x_update,J_update,filter.track_threshold);  T_prune= length(r_update);
    [r_update,w_update,x_update,J_update]= cap_tracks(r_update,w_update,x_update,J_update,filter.T_max);              T_cap  = length(r_update);  
    
    T_update= T_cap;
    
    %(within tracks) resampling
    for t=1:T_update
        J_rsp(t)= max(round(r_update(t)*filter.J_max),filter.J_min);
        idx_rsp= randsample(length(w_update{t}),J_rsp(t),true,w_update{t}); %idx_rsp= resample(w_update{t},J_rsp(t));
        w_update{t}= ones(J_rsp(t),1)/J_rsp(t);
        x_update{t}= x_update{t}(:,idx_rsp);
    end
    
    J_update= J_rsp;
    
    %--- state extraction
    cdn_update= prod(1-r_update)*esf(r_update./(1-r_update));
    [~,idx_max_cdn] = max(cdn_update);
    map_cdn = idx_max_cdn-1;
    est.N(k) = min(length(r_update),map_cdn);
    [~,idx_max_tgt]= sort(-r_update);
    for t=1:est.N(k)
        est.X{k} = [est.X{k} x_update{idx_max_tgt(t)}*w_update{idx_max_tgt(t)}];
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
         ' #samp filt=',num2str(sum(J_rsp))   ]);
    end

end

function clipped_r= limit_range(r)

r(r>0.999)=0.999;
r(r<0.001)=0.001;
clipped_r= r;

function [new_r,new_w,new_x,new_J]= prune_tracks(r,w,x,J,track_threshold)

idx= find(r>track_threshold);

new_r= zeros(length(idx));
new_w= cell(length(idx),1);
new_x= cell(length(idx),1);
new_J= zeros(length(idx),1);

new_r= r(idx);
for i=1:length(idx)
    new_w{i}= w{idx(i)};
    new_x{i}= x{idx(i)};
end
new_J= J(idx);

function [new_r,new_w,new_x,new_J]= cap_tracks(r,w,x,J,T_max)

if length(r) > T_max

[~,idx]= sort(-r);

new_r= zeros(length(idx));
new_w= cell(length(idx),1);
new_x= cell(length(idx),1);
new_J= zeros(length(idx),1);

new_r= r(idx(1:T_max));
for i=1:T_max
    new_w{i}= w{idx(i)};
    new_x{i}= m{idx(i)};
end
new_J= J(idx(1:T_max));

else
    new_r= r;
    new_w= w;
    new_x= x;
    new_J= J;
end


            