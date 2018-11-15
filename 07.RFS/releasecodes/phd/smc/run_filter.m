function est = run_filter(model,meas)

% This is the MATLAB code for the SMC-PHD filter proposed in
% (assuming no target spawning)
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
filter.J_target= 1000;                                        %generated number of particles per expected target
filter.J_birth= model.L_birth*filter.J_target;                %generated number of particles from birth intensity

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
w_update= eps;
m_init= [0.1;0;0.1;0;0.01];
P_init= diag([100 10 100 10 1]).^2;
x_update= m_init;

%recursive filtering
for k=1:meas.K
    %---prediction 
    pS_vals= compute_pS(model,x_update); pS_vals= pS_vals(:);
    qS_vals= 1-pS_vals;

    x_predict = gen_newstate_fn(model,x_update,'noise');                                                                %surviving particles
    w_predict = pS_vals.*w_update;                                                                                      %surviving weights

    x_predict= cat(2,gen_gms(model.w_birth,model.m_birth,model.P_birth,filter.J_birth),x_predict);                      %append birth particles
    w_predict= cat(1,sum(model.w_birth)*ones(filter.J_birth,1)/filter.J_birth,w_predict);                               %append birth weights                                                
        
    %---update
    %number of measurements
    m= size(meas.Z{k},2);
    pD_vals= compute_pD(model,x_predict); pD_vals= pD_vals(:);
    qD_vals= 1-pD_vals;

    %missed detection weight
    pseudo_likelihood= qD_vals;
    
    if m~=0
        %m detection weights
        for ell=1:m
            meas_likelihood= compute_likelihood(model,meas.Z{k}(:,ell),x_predict)';
            pseudo_likelihood = pseudo_likelihood+pD_vals.*meas_likelihood/(model.lambda_c*model.pdf_c + (pD_vals.*meas_likelihood)'*w_predict);
        end
    end
    w_update= pseudo_likelihood.*w_predict;
    x_update= x_predict;    
            
    %---for diagnostics
    w_posterior= w_update;
    
    %---resampling
    J_rsp= min(ceil(sum(w_update)*filter.J_target),filter.J_max);
    idx= randsample(length(w_update),J_rsp,true,w_update); %idx= resample(w_update,J_rsp);
    w_update= sum(w_update)*ones(J_rsp,1)/J_rsp;
    x_update= x_update(:,idx);
 
    %--- state extraction   
    if sum(w_update) > 0.5
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
  
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #est mean=' num2str(sum(w_update),4),...
         ' #est card=' num2str(est.N(k),4),...
         ' Neff_updt= ',num2str(round(1/sum((w_posterior/sum(w_posterior)).^2)))...
         ' Neff_rsmp= ',num2str(round(1/sum((w_update/sum(w_update)).^2)))   ]);
    end

end

            