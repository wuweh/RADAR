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
filter.L_max= 100;                  %limit on number of Gaussians
filter.elim_threshold= 1e-5;        %pruning threshold
filter.merge_threshold= 4;          %merging threshold

filter.P_G= 0.9999999;                           %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value
filter.gate_flag= 1;                             %gating on or off 1/0

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
w_update(1)= 1;
m_update(:,1)= [0.1;0;0.1;0;0.01];
P_update(:,:,1)= diag([100 10 100 10 1]).^2;
L_update = 1;

%recursive filtering
for k=1:meas.K
    %---prediction 
    [m_predict,P_predict] = ekf_predict_multiple(model,m_update,P_update);
    w_predict= w_update;
    L_predict= L_update;
    
    %---gating
    if filter.gate_flag
        meas.Z{k}= gate_meas_ekf(meas.Z{k},filter.gamma,model,m_predict,P_predict);     
    end
        
    %---update
    %number of measurements
    m= size(meas.Z{k},2);
    
    %missed detection term - scale to get original expression with factor exp(-model.lambda_c)*(model.lambda_c)^(m-1)*(model.pdf_c)^(m-1)
    w_update = model.Q_D*w_predict*(model.lambda_c)*(model.pdf_c);
    m_update = m_predict;
    P_update = P_predict;
    
    if m~=0
        %m detection terms - scale to get original expression with factor exp(-model.lambda_c)*(model.lambda_c)^(m-1)*(model.pdf_c)^(m-1)
        [qz_temp,m_temp,P_temp] = ekf_update_multiple(meas.Z{k},model,m_predict,P_predict);
        for ell=1:m
            w_temp = model.P_D*w_predict(:).*qz_temp(:,ell);
            w_update = cat(1,w_update,w_temp);
            m_update = cat(2,m_update,m_temp(:,:,ell));
            P_update = cat(3,P_update,P_temp);
        end
    end
            
    %normalize weights
    w_update = w_update/sum(w_update);        
            
    %---mixture management
    L_posterior= length(w_update);
    
    %pruning, merging, capping
    [w_update,m_update,P_update]= gaus_prune(w_update,m_update,P_update,filter.elim_threshold);    L_prune= length(w_update);
    [w_update,m_update,P_update]= gaus_merge(w_update,m_update,P_update,filter.merge_threshold);   L_merge= length(w_update);
    [w_update,m_update,P_update]= gaus_cap(w_update,m_update,P_update,filter.L_max);               L_cap  = length(w_update);
    
    L_update= L_cap;
    
    %--- state extraction
    [~,idx]= max(w_update);
    est.X{k} = m_update(:,idx);
    est.N(k)= 1;
    est.L{k}= [];
  
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #gaus orig=',num2str(L_posterior),...
         ' #gaus elim=',num2str(L_prune), ...
         ' #gaus merg=',num2str(L_merge)   ]);
    end

end

            