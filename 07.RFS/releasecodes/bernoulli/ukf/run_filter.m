function est = run_filter(model,meas)

% This is the MATLAB code for the Bernoulli filter with RFS observations proposed in
% (for a single sensor only)
% B.-T. Vo, C.M. See, N. Ma and W.T. Ng, "Multi-Sensor Joint Detection and Tracking with the Bernoulli Filter," IEEE Trans. Aerospace and Electronic Systems, Vol. 48, No. 2, pp. 1385 - 1402, 2012.
% http://ba-ngu.vo-au.com/vo/VSMN_Bernoulli_TAES12.pdf
% ---BibTeX entry
% @ARTICLE{BER,
% author={B.-T. Vo and C.M. See and N. Ma and W.T. Ng},
% journal={IEEE Transactions on Aerospace and Electronic Systems},
% title={Multi-Sensor Joint Detection and Tracking with the Bernoulli Filter},
% year={2012},
% month={April},
% volume={48},
% number={2},
% pages={1385-1402}}
%---
% ... see also ...
%---
% B. Ristic, B.-T. Vo, B.-N. Vo, and A. Farina "A Tutorial on Bernoulli Filters: Theory, Implementation and Applications," IEEE Trans. Signal Processing, Vol. 61, No. 13, pp. 3406 - 3430, 2013.
% http://ba-ngu.vo-au.com/vo/RVVF_Bernoulli_TSP13.pdf
% ---BibTeX entry
% @ARTICLE{BERTUT,
% author={B. Ristic and B.-T. Vo and B.-N. Vo and A. Farina},
% journal={IEEE Transactions on Signal Processing},
% title={A Tutorial on Bernoulli Filters: Theory, Implementation and Applications},
% year={2013},
% month={July},
% volume={61},
% number={13},
% pages={3406-3430}}
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
filter.gate_flag= 0;                             %gating on or off 1/0

filter.ukf_alpha= 1;                %scale parameter for UKF - choose alpha=1 ensuring lambda=beta and offset of first cov weight is beta for numerical stability
filter.ukf_beta= 2;                 %scale parameter for UKF
filter.ukf_kappa= 2;                %scale parameter for UKF (alpha=1 preferred for stability, giving lambda=1, offset of beta for first cov weight)
filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
r_update= 0.001;
w_update(1)= 1;
m_update(:,1)= [0.1;0;0.1;0;0.01];
P_update(:,:,1)= diag([100 10 100 10 1]).^2;
L_update = 1;

%recursive filtering
for k=1:meas.K
    %---prediction 
    r_predict= model.r_birth*(1-r_update) + model.P_S*r_update;                                     %predicted existence probability

    [m_predict,P_predict] = ukf_predict_multiple(model,m_update,P_update,filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);                          %surviving components
    w_predict= model.P_S*r_update*w_update;                                                                                                           %surviving weights

    m_predict= cat(2,model.m_birth{1},m_predict); P_predict=cat(3,model.P_birth{1},P_predict);      %append birth components
    w_predict= cat(1,model.r_birth*(1-r_update)*model.w_birth{1},w_predict);                        %append birth weights
                                                 
    w_predict= w_predict/sum(w_predict);                                                            %normalize weights
    r_predict= limit_range(r_predict);                                                              %limit range of 0<r<1 for numerical stability
    L_predict= model.L_birth+L_update;                                                              %number of predicted components
    
    %---gating
    if filter.gate_flag
        meas.Z{k}= gate_meas_ukf(meas.Z{k},filter.gamma,model,m_predict,P_predict,filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);       
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
        [qz_temp,m_temp,P_temp] = ukf_update_multiple(meas.Z{k},model,m_predict,P_predict,filter.ukf_alpha,filter.ukf_kappa,filter.ukf_beta);
        for ell=1:m
            w_temp = model.P_D*w_predict(:).*qz_temp(:,ell);
            w_update = cat(1,w_update,w_temp);
            m_update = cat(2,m_update,m_temp(:,:,ell));
            P_update = cat(3,P_update,P_temp);
        end
    end

    %existence probability
    r_update= (r_predict*sum(w_update))/( (model.lambda_c*model.pdf_c)*(1-r_predict) + r_predict*sum(w_update) );
    r_update= limit_range(r_update); 
            
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
    if r_update>0.5
        [~,idx]= max(w_update);
        est.X{k} = m_update(:,idx);
        est.N(k)= 1;
        est.L{k}= [];
    end
  
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #r prob=' num2str(r_update,4),...
         ' #est card=' num2str(est.N(k),4),...
         ' #gaus orig=',num2str(L_posterior),...
         ' #gaus elim=',num2str(L_prune), ...
         ' #gaus merg=',num2str(L_merge)   ]);
    end

end

function clipped_r= limit_range(r)

r(r>0.999)=0.999;
r(r<0.001)=0.001;
clipped_r= r;


            