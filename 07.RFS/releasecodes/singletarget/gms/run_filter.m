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

filter.P_G= 0.999;                               %gate size in percentage
filter.gamma= chi2inv(filter.P_G,model.z_dim);   %inv chi^2 dn gamma value  卡方门限值
filter.gate_flag= 1;                             %gating on or off 1/0

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering 

%initial prior
w_update(1)= 1;
m_update(:,1)= [0.1;0;0.1;0];
P_update(:,:,1)= diag([1000 100 1000 100]).^2;
L_update = 1;

%recursive filtering
for k=1:meas.K
    %---prediction 
    %得到当前目标量测预测和协方差预测
    [m_predict,P_predict] = kalman_predict_multiple(model,m_update,P_update);
    w_predict= w_update;
    L_predict= L_update;
    
    %---gating
    %根据卡方门限值来选择进入门限的两侧目标
    if filter.gate_flag
        meas.Z{k}= gate_meas_gms(meas.Z{k},filter.gamma,model,m_predict,P_predict);        
    end
        
    %---update
    %number of measurements
    %获得当前落入门限的量测目标个数
    m= size(meas.Z{k},2);
    
    %missed detection term - scale to get original expression with factor exp(-model.lambda_c)*(model.lambda_c)^(m-1)*(model.pdf_c)^(m-1)
    %===================对漏检目标更新======================
    % model.Q_D : 漏检概率
    % w_predict 上一时刻的预测强度
    % lambda_c  poisson average rate of uniform clutter (per scan)
    % pdf_c 杂波密度
    
    %无目标时的 w m p
    w_update = model.Q_D*w_predict*(model.lambda_c)*(model.pdf_c);
    m_update = m_predict;
    P_update = P_predict;
    
    %===================对真实目标更新======================
    if m~=0
        %m detection terms - scale to get original expression with factor exp(-model.lambda_c)*(model.lambda_c)^(m-1)*(model.pdf_c)^(m-1)
        [qz_temp,m_temp,P_temp] = kalman_update_multiple(meas.Z{k},model,m_predict,P_predict);
        for ell=1:m
            w_temp = model.P_D*w_predict(:).*qz_temp(:,ell);
            w_update = cat(1,w_update,w_temp);
            m_update = cat(2,m_update,m_temp(:,:,ell));         %后验状态集
            P_update = cat(3,P_update,P_temp);                      %后验协方差集
        end
    end
            
    %normalize weights
    w_update = w_update/sum(w_update);        
            
    %---mixture management
    L_posterior= length(w_update);
    
    %pruning, merging, capping
    %高斯修剪融合算法
    [w_update,m_update,P_update]= gaus_prune(w_update,m_update,P_update,filter.elim_threshold);    L_prune= length(w_update);                   %只选取大于门限值的目标
    [w_update,m_update,P_update]= gaus_merge(w_update,m_update,P_update,filter.merge_threshold);   L_merge= length(w_update);               % 融合门限merge_threshold
    [w_update,m_update,P_update]= gaus_cap(w_update,m_update,P_update,filter.L_max);               L_cap  = length(w_update);
    
    L_update= L_cap;
    
    %--- state extraction
    [~,idx]= max(w_update);  %选择权重最大的高斯项
    est.X{k} = m_update(:,idx); %更新预测值
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

            