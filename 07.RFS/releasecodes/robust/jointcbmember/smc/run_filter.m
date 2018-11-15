function est = run_filter(model,meas)

% This is the MATLAB code for the Robust CBMeMBer filter proposed in
% (without track labelling)
% B.-T. Vo, B.-N. Vo, R. Hoseinnezhad, and R. Mahler "Robust Multi-Bernoulli Filtering," IEEE Journal on Selected Topics in Signal Processing, Vol. 7, No. 3, pp. 399-409, 2013.
% http://ba-ngu.vo-au.com/vo/VVHM_JSSP13.pdf
% ---BibTeX entry
% @ARTICLE{RobustCBMEMBER,
% author={B.-T. Vo and B.-N. Vo and R. Hoseinnezhad and R. Mahler},
% journal={IEEE Transactions on Signal Processing},
% title={Robust Multi-Bernoulli Filtering},
% year={2013},
% month={Jun},
% volume={7},
% number={3},
% pages={399-409}} 
%---

%=== Setup

%output variables
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);
est.L= zeros(meas.K,1);
est.pD= cell(meas.K,1);

%filter parameters
filter.T_max= 100;                  %maximum number of tracks
filter.track_threshold= 1e-3;       %threshold to prune tracks

filter.J_max= 1000;                  %max number of particles in each track (otherwise allocated proportion to existence probability)
filter.J_min= 300;                   %min number of particles in each track (otherwise allocated proportion to existence probability)

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior

T_update= 10*size(meas.Z(1),2);
r_update= 1/T_update*ones(T_update,1);
J_update= max(round(r_update*filter.J_max/2),round(filter.J_min/2));
x_init= [0.1;0;0.1;0;0.01];
P_init= diag([100 10 100 10 1]).^2;
for j=1:T_update
   w_update{j}= ones(J_update(j),1)*1/J_update(j);
   x_update{j}= [ betarnd(5,5,1,J_update(j)); repmat(x_init, [1, J_update(j)])+chol(P_init)*randn(model.x_dim-1,J_update(j))];
   l_update{j}= zeros(J_update(j),1);  
end


%recursive filtering
for k=1:meas.K
    %---prediction 
    clear r_predict J_predict w_predict x_predict l_predict;
    T_predict= T_update+model.T_birth+model.T_birth_clt;                                                                              %total number of tracks/components
    J_predict= zeros(T_predict,1);                                                                                  %total number of particles in each track/component
    r_predict= zeros(T_predict,1);                                                                                  %existence probability for tracks/components
    x_predict= cell(T_predict,1);                                                                                   %states of particles in each track/component
    w_predict= cell(T_predict,1);                                                                                   %weights of particles in each track/component
    l_predict= cell(T_predict,1);
    
    r_predict= model.r_birth;
    for t=1:model.T_birth
        J_predict(t)= max(round(model.r_birth(t)*filter.J_max),filter.J_min);                                                                            %append birth particles count
%         x_birth_temp1= gen_gms(model.w_birth{t},model.m_birth{t},model.P_birth{t},J_predict(t));
        x_birth_temp1= repmat(model.m_birth{t}, [1, J_predict(t)])+model.B_birth{t}*randn(model.x_dim-1,J_predict(t));
        x_predict{t}= [betarnd(model.u_b_tg(t),model.v_b_tg(t),1,J_predict(t)); x_birth_temp1];                                                   %append birth particles
        w_predict{t}= ones(J_predict(t),1)/J_predict(t);                                                                                   %append birth weights
        l_predict{t}= ones(J_predict(t),1);
    end     
    
    offset= model.T_birth;
    r_predict= [r_predict; sum(model.r_birth_clt)/model.T_birth_clt*ones(model.T_birth_clt,1)];
%     r_predict= [r_predict; model.r_birth_clt];
    for t=1:model.T_birth_clt             
        J_predict(offset+t)= max(round(model.r_birth_clt(t)*filter.J_max/2),round(filter.J_min/2));                                                                            %append birth particles count
        x_birth_temp0= repmat(model.m_birth_clt{t}, [1, J_predict(offset+t)])+model.B_birth_clt{t}*randn(model.x_dim-1,J_predict(offset+t));
        x_predict{offset+t}= [betarnd(model.u_b_clt(t),model.v_b_clt(t),1,J_predict(offset+t)); x_birth_temp0];                                                   %append birth particles
        w_predict{offset+t}= ones(J_predict(offset+t),1)/J_predict(offset+t);
        l_predict{offset+t}= zeros(J_predict(offset+t),1);
    end
    
    offset= model.T_birth+model.T_birth_clt;
    %prediction for surviving tracks
    for t=1:T_update
        idx1= find(l_update{t} == 1);
        if ~isempty(idx1)
            w_temp1= w_update{t}(idx1).*compute_pS_tg(model,x_update{t}(:,idx1));
        else
            w_temp1=[];
        end
        idx0= find(l_update{t} == 0);
        if ~isempty(idx0)
            w_temp0= w_update{t}(idx0).*compute_pS_clt(model,x_update{t}(:,idx0));
        else
            w_temp0= [];
        end
        
        r_predict(offset+t)= r_update(t)* (sum(w_temp1)+sum(w_temp0));
        J_predict(offset+t)= J_update(t);
        w_predict{offset+t}= [w_temp1/(sum(w_temp1)+sum(w_temp0)); w_temp0/(sum(w_temp1)+sum(w_temp0))];
        x_predict{offset+t}= [gen_newstate_tg(model,x_update{t}(:,idx1)) gen_newstate_clt(model,x_update{t}(:,idx0))];
        l_predict{offset+t}= [ones(length(idx1),1); zeros(length(idx0),1) ];
    end

    r_predict= limit_range(r_predict);                                                                         %limit range of 0<r<1 for numerical stability

    
    
    %---construction of pseudo PHD for update
    J_pseudo= sum(J_predict);                                           %number of particles in pseudo-PHD
    x_pseudo= zeros(model.x_dim,J_pseudo);                              %states of particles in pseudo-PHD
    l_pseudo= zeros(J_pseudo,1);
    w_pseudo= zeros(J_pseudo,1);                                        %weights of particles in pseudo-PHD
    w_pseudo1= zeros(J_pseudo,1);                                       %alt weight (1) of particles in pseudo-PHD - used in CB-MeMBer update later
    w_pseudo2= zeros(J_pseudo,1);                                       %alt weight (2) of particles in pseudo-PHD - used in CB-MeMBer update later
    
    start_pt= 1;
    for t=1:T_predict                  
        end_pt= start_pt+J_predict(t)-1;
        idx1= find(l_predict{t}==1);
        if ~isempty(idx1)
            w_temp1= w_predict{t}(idx1).*x_predict{t}(1,idx1)';
            x_temp1= x_predict{t}(:,idx1);
        else
            w_temp1= [];
            x_temp1= [];
        end
        idx0= find(l_predict{t}==0);
        if ~isempty(idx0)
            w_temp0= w_predict{t}(idx0).*x_predict{t}(1,idx0)';
            x_temp0= x_predict{t}(:,idx0);
        else
            w_temp0= [];
            x_temp0= [];
        end
        
        x_pseudo(:,start_pt:end_pt)= cat(2,x_temp1, x_temp0);
        l_pseudo(start_pt:end_pt)= cat(1,ones(length(idx1),1), zeros(length(idx0),1));
        w_pseudo(start_pt:end_pt) = cat(1,r_predict(t)/(1-r_predict(t))*w_predict{t}(idx1), r_predict(t)/(1-r_predict(t))*w_predict{t}(idx0));
        w_pseudo1(start_pt:end_pt)= cat(1,r_predict(t)/(1-r_predict(t)*(sum(w_temp1)+sum(w_temp0)))*w_predict{t}(idx1), r_predict(t)/(1-r_predict(t)*(sum(w_temp1)+sum(w_temp0)))*w_predict{t}(idx0));
        w_pseudo2(start_pt:end_pt)= cat(1,r_predict(t)*(1-r_predict(t))/((1-r_predict(t)*(sum(w_temp1)+sum(w_temp0)))^2)*w_predict{t}(idx1), r_predict(t)*(1-r_predict(t))/((1-r_predict(t)*(sum(w_temp1)+sum(w_temp0)))^2)*w_predict{t}(idx0));

        start_pt= end_pt+1;
    end
    %reorder to adhere to tag convention - tag 1 particles, tag 0 particles
    idx1= find(l_pseudo==1); 
    idx0= find(l_pseudo==0); 
    x_pseudo= [x_pseudo(:,idx1) x_pseudo(:,idx0)];
    l_pseudo= [l_pseudo(idx1);l_pseudo(idx0)];
    w_pseudo= [w_pseudo(idx1);w_pseudo(idx0)];
    w_pseudo1= [w_pseudo1(idx1);w_pseudo1(idx0)];
    w_pseudo2= [w_pseudo2(idx1);w_pseudo2(idx0)];
    
    
    %---update
    %number of measurements
    m= size(meas.Z{k},2);
    clear r_update J_update w_update x_update l_update r_comp1 r_comp0;
    
    T_update= T_predict+m;                                                                                        %total number of tracks/components
    J_update= zeros(T_update,1);                                                                                  %total number of particles in each track/component
    r_update= zeros(T_update,1);                                                                                  %existence probability for tracks/components
    x_update= cell(T_update,1);                                                                                   %states of particles in each track/component
    w_update= cell(T_update,1);                                                                                   %weights of particles in each track/component
 
    %legacy tracks
    for t=1:T_predict  
        idx1= find(l_predict{t}==1);
        if ~isempty(idx1)
            w_temp1= w_predict{t}(idx1).*x_predict{t}(1,idx1)';
            m_temp1= w_predict{t}(idx1).*(1-x_predict{t}(1,idx1)');
        else
            w_temp1= [];
            m_temp1= [];
        end   
        
        idx0= find(l_predict{t}==0);
        if ~isempty(idx0)
            w_temp0= w_predict{t}(idx0).*x_predict{t}(1,idx0)';
            m_temp0= w_predict{t}(idx0).*(1-x_predict{t}(1,idx0)');
        else
            w_temp0= [];
            m_temp0= [];
        end
        
        J_update(t)= J_predict(t);
        rtrack1= r_predict(t)*sum(m_temp1)/(1-r_predict(t)*(sum(w_temp1)+sum(w_temp0)));
        rtrack0= r_predict(t)*sum(m_temp0)/(1-r_predict(t)*(sum(w_temp1)+sum(w_temp0)));
        r_update(t)= rtrack1+rtrack0;
        r_comp1(t)= rtrack1;
        r_comp0(t)= rtrack0;
        
        w_update{t}= [m_temp1/(sum(m_temp1)+sum(m_temp0)) ; m_temp0/(sum(m_temp1)+sum(m_temp0))]; 
        x_update{t}= x_predict{t};
        l_update{t}= l_predict{t};
        
    end
    %measurement updated tracks
    if m~=0
        offset= T_predict;
        idx1= find(l_pseudo==1);
        idx0= find(l_pseudo==0);

        pD_vals_tg= x_pseudo(1,idx1)';
        pD_vals_clt= x_pseudo(1,idx0)';
        for ell=1:m

            meas_likelihood_tg= compute_likelihood_tg(model,meas.Z{k}(:,ell),x_pseudo(:,idx1))';
            meas_likelihood_clt= compute_likelihood_clt(model,meas.Z{k}(:,ell),x_pseudo(:,idx0))';
          
            r_temp1= sum(pD_vals_tg.*meas_likelihood_tg.*w_pseudo2(idx1))/(sum(pD_vals_tg.*meas_likelihood_tg.*w_pseudo1(idx1))+sum(pD_vals_clt.*meas_likelihood_clt.*w_pseudo1(idx0)));
            r_temp0= sum(pD_vals_clt.*meas_likelihood_clt.*w_pseudo2(idx0))/(sum(pD_vals_tg.*meas_likelihood_tg.*w_pseudo1(idx1))+sum(pD_vals_clt.*meas_likelihood_clt.*w_pseudo1(idx0)));
            r_track= r_temp1+r_temp0;
            
            w_temp1= pD_vals_tg.*meas_likelihood_tg.*w_pseudo(idx1);
            w_temp0= pD_vals_clt.*meas_likelihood_clt.*w_pseudo(idx0);
                     
            r_update(offset+ell)= r_track;
            w_update{offset+ell}= [w_temp1/(sum(w_temp1)+sum(w_temp0)); w_temp0/(sum(w_temp1)+sum(w_temp0))];
            J_update(offset+ell)= J_pseudo;
            x_update{offset+ell}= x_pseudo;
            l_update{offset+ell}= l_pseudo;
            
            r_comp1(offset+ell)= r_temp1;
            r_comp0(offset+ell)= r_temp0;
    
        end
        
    end
    r_update= limit_range(r_update);
            
    %---track/component management
    T_posterior= T_update;
    
   
    %(within tracks) resampling
    for t=1:T_update
        J_rsp(t)= max(round(r_update(t)*filter.J_max),filter.J_min);
        idx_rsp= randsample(length(w_update{t}),J_rsp(t),true,w_update{t}); %idx_rsp= resample(w_update{t},J_rsp(t));
        w_update{t}= ones(J_rsp(t),1)/J_rsp(t);
        x_update{t}= x_update{t}(:,idx_rsp);
        l_update{t}= l_update{t}(idx_rsp);
    end
    J_update= J_rsp;
        
    %tracks: pruning, capping
    [r_update,r_comp1,r_comp0,w_update,x_update,J_update,l_update]= prune_tracks(r_update,r_comp1,r_comp0,w_update,x_update,J_update,l_update,filter.track_threshold);  T_prune= length(r_update);
    [r_update,r_comp1,r_comp0,w_update,x_update,J_update,l_update]= cap_tracks(r_update,r_comp1,r_comp0,w_update,x_update,J_update,l_update,filter.T_max);              T_cap  = length(r_update);  
    
    T_update= T_cap;

        
    %--- state extraction
    pD_tmp= [];
    est.soft_N(k)= sum(r_comp1);
    est.N(k)= round(est.soft_N(k));
    [~,idx]= sort(-r_comp1);
    for t=1:min(est.N(k),T_update)
        idx1= find(l_update{idx(t)}==1);
        est.X{k} = [est.X{k} x_update{idx(t)}(:,idx1)*w_update{idx(t)}(idx1)];
        pD_tmp = [pD_tmp; x_update{idx(t)}(1,idx1)*w_update{idx(t)}(idx1)];
    end
    est.pD{k} = mean(pD_tmp);   
    
    %estimate clutter rate
    avg_pd0= zeros(T_update,1);
    for j=1:T_update
        idx0= find(l_update{j}==0);
        if ~isempty(idx0)
            avg_pd0(j)= (w_update{j}(idx0))'*(x_update{j}(1,idx0))';
        end
    end;
    est.L(k)= sum(r_comp0'.*avg_pd0);   
  
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #avg pD=' num2str(est.pD{k}),...         
         ' #avg lambda=' num2str(est.L(k)),...
         ' #eap target=' num2str(sum(r_comp1)),...
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

function [new_r,new_r_com1,new_r_com0,new_w,new_x,new_J,new_l]= prune_tracks(r,r_com1,r_com0,w,x,J,l,track_threshold)

idx= find(r>track_threshold);

new_r= zeros(length(idx));
new_r_com1= zeros(length(idx));
new_r_com0= zeros(length(idx));
new_w= cell(length(idx),1);
new_x= cell(length(idx),1);
new_J= zeros(length(idx),1);
new_l= cell(length(idx),1);

new_r= r(idx); new_r_com1= r_com1(idx); new_r_com0= r_com0(idx);
for i=1:length(idx)
    new_w{i}= w{idx(i)};
    new_x{i}= x{idx(i)};
    new_l{i}= l{idx(i)};
end
new_J= J(idx);

function [new_r,new_r_com1,new_r_com0,new_w,new_x,new_J,new_l]= cap_tracks(r,r_com1,r_com0,w,x,J,l,T_max)

if length(r) > T_max

[~,idx]= sort(-r);

new_r= zeros(length(idx));
new_r_com1= zeros(length(idx));
new_r_com0= zeros(length(idx));
new_w= cell(length(idx),1);
new_x= cell(length(idx),1);
new_J= zeros(length(idx),1);
new_l= cell(length(idx),1);

new_r= r(idx(1:T_max)); new_r_com1= r_com1(idx(1:T_max)); new_r_com0= r_com0(idx(1:T_max));
for i=1:T_max
    new_w{i}= w{idx(i)};
    new_x{i}= x{idx(i)};
    new_l{i}= l{idx(i)};
end
new_J= J(idx(1:T_max));

else
    new_r= r;
    new_r_com1= r_com1;
    new_r_com0= r_com0;
    new_w= w;
    new_x= x;
    new_J= J;
    new_l= l;
end


            