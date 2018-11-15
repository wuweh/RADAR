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
filter.J_max= 3000;                  %total number of particles
filter.J_birth= 1000;                %total number of birth particles only

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
r_update= 0.001;
w_update= ones(filter.J_max,1)/filter.J_max;
m_init= [0.1;0;0.1;0;0.01];
P_init= diag([100 10 100 10 1]).^2;
x_update= gen_gms(1,m_init,P_init,filter.J_max);

%recursive filtering
for k=1:meas.K
    %---prediction 
    pS_vals= compute_pS(model,x_update); pS_vals= pS_vals(:);
    qS_vals= 1-pS_vals;
    
    r_predict= model.r_birth*(1-r_update) + pS_vals'*w_update*r_update;                                                 %predicted existence probability

    x_predict = gen_newstate_fn(model,x_update,'noise');                                                                %surviving particles
    w_predict = r_update*pS_vals.*w_update;                                                                             %surviving weights

    x_predict= cat(2,gen_gms(model.w_birth{1},model.m_birth{1},model.P_birth{1},filter.J_birth),x_predict);             %append birth components
    w_predict= cat(1,model.r_birth*(1-r_update)*ones(filter.J_birth,1)/filter.J_birth,w_predict);                       %append birth weights
                                                 
    w_predict= w_predict/sum(w_predict);                                                                                %normalize weights
    r_predict= limit_range(r_predict);                                                                                  %limit range of 0<r<1 for numerical stability
        
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

    %existence probability
    r_update= (r_predict*sum(w_update))/( (model.lambda_c*model.pdf_c)*(1-r_predict) + r_predict*sum(w_update) );
    r_update= limit_range(r_update);
            
    %normalize weights
    w_update = w_update/sum(w_update);        
            
    %---for diagnostics
    w_posterior= w_update;
    
    %---resampling
    idx= randsample(length(w_update),filter.J_max,true,w_update); %idx= resample(w_update,filter.J_max);
    w_update= ones(filter.J_max,1)/filter.J_max;
    x_update= x_update(:,idx);
    
    %--- state extraction
    if r_update>0.5
        est.X{k} = x_update*w_update;
        est.N(k)= 1;
        est.L{k}= [];
    end
  
    %---display diagnostics
    if ~strcmp(filter.run_flag,'silence')
        disp([' time= ',num2str(k),...
         ' #r prob=' num2str(r_update,4),...
         ' #est card=' num2str(est.N(k),4),...
         ' Neff_updt= ',num2str(round(1/sum(w_posterior.^2)))...
         ' Neff_rsmp= ',num2str(round(1/sum(w_update.^2)))   ]);
    end

end

function clipped_r= limit_range(r)

r(r>0.999)=0.999;
r(r<0.001)=0.001;
clipped_r= r;


            