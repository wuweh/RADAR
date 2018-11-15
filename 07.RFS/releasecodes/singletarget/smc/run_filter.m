function est = run_filter(model,meas)

% This is the MATLAB code for the single target filter with RFS observations proposed in
% (assuming Poisson clutter and no extraneous measurements)
% B.-T. Vo, B.-N. Vo, and A. Cantoni, "Bayesian filtering with random finite set observations," IEEE Trans. Signal Processing, Vol. 56, No. 4, pp. 1313-1326, 2008.
% http://ba-ngu.vo-au.com/vo/VVCsingletargetSP08.pdf
% ---BibTeX entry
% @ARTICLE{STF, 
% author={B.-T.Vo and B.-N. Vo and A. Cantoni},
% journal={IEEE Transactions on Signal Processing},
% title={Bayesian Filtering With Random Finite Set Observations},
% year={2008},
% month={April},
% volume={56},
% number={4},
% pages={1313-1326}}
%---

%=== Setup

%output variables
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);
est.L= cell(meas.K,1);

%filter parameters
filter.J_max= 3000;                  %total number of particles

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
w_update= ones(filter.J_max,1)/filter.J_max;
m_init= [0.1;0;0.1;0;0.01];
P_init= diag([100 10 100 10 1]).^2;
x_update= gen_gms(1,m_init,P_init,filter.J_max);

%recursive filtering
for k=1:meas.K
    %---prediction 
    x_predict = gen_newstate_fn(model,x_update,'noise');
    w_predict= w_update;
        
    %---update
    %number of measurements
    m= size(meas.Z{k},2);
    pD_vals= compute_pD(model,x_predict); pD_vals= pD_vals(:);
    qD_vals= 1-pD_vals;

    %missed detection weight - scale to get original expression with factor exp(-model.lambda_c)*(model.lambda_c)^(m-1)*(model.pdf_c)^(m-1) 
    rfs_likelihood= qD_vals*(model.lambda_c)*(model.pdf_c);
    
    if m~=0
        %m detection weights - scale to get original expression with factor exp(-model.lambda_c)*(model.lambda_c)^(m-1)*(model.pdf_c)^(m-1)
        for ell=1:m
            rfs_likelihood = rfs_likelihood+pD_vals.*compute_likelihood(model,meas.Z{k}(:,ell),x_predict)';
        end
    end
    w_update= rfs_likelihood.*w_predict;
    x_update= x_predict;
            
    %normalize weights
    w_update = w_update/sum(w_update);        
            
    %---for diagnostics
    w_posterior= w_update;
    
    %---resampling
    idx= randsample(length(w_update),filter.J_max,true,w_update); %idx= resample(w_update,filter.J_max);
    w_update= ones(filter.J_max,1)/filter.J_max;
    x_update= x_update(:,idx);
    
    %--- state extraction
    est.X{k} = x_update*w_update;
    est.N(k)= 1;
    est.L{k}= [];
  
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' Neff_updt= ',num2str(round(1/sum(w_posterior.^2)))...
         ' Neff_rsmp= ',num2str(round(1/sum(w_update.^2)))   ]);
    end

end

            