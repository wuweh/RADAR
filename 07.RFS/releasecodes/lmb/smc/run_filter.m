function est = run_filter(model,meas)

% This is the MATLAB code for the Labeled Multi-Bernoulli filter proposed in
% S. Reuter, B.-T. Vo, B.-N. Vo, and K. Dietmayer, "The labelled multi-Bernoulli filter," IEEE Trans. Signal Processing, Vol. 62, No. 12, pp. 3246-3260, 2014
% http://ba-ngu.vo-au.com/vo/RVVD_LMB_TSP14.pdf
% which propagates an LMB approximation of the GLMB update proposed in
% B.-N. Vo, B.-T. Vo, and D. Phung, "Labeled Random Finite Sets and the Bayes Multi-Target Tracking Filter," IEEE Trans. Signal Processing, Vol. 62, No. 24, pp. 6554-6567, 2014
% http://ba-ngu.vo-au.com/vo/VVP_GLMB_TSP14.pdf
%
% Note 1: no dynamic grouping or adaptive birth is implemented in this code, only the standard filter with static birth is given
% Note 2: the simple example used here is the same as in the CB-MeMBer filter code for a quick demonstration and comparison purposes
% Note 3: more difficult scenarios require more components/hypotheses (thus exec time) and/or a better lookahead
% ---BibTeX entry
% @ARTICLE{LMB,
% author={S. Reuter and B.-T. Vo and B.-N. Vo and K. Dietmayer},
% journal={IEEE Transactions on Signal Processing},
% title={The Labeled Multi-Bernoulli Filter},
% year={2014},
% month={Jun}
% volume={62},
% number={12},
% pages={3246-3260}}
%
% @ARTICLE{GLMB2,
% author={B.-T. Vo and B.-N. Vo and D. Phung},
% journal={IEEE Transactions on Signal Processing},
% title={Labeled Random Finite Sets and the Bayes Multi-Target Tracking Filter},
% year={2014},
% month={Dec}
% volume={62},
% number={24},
% pages={6554-6567}}
%---

%=== Setup

%output variables
est.X= cell(meas.K,1);
est.N= zeros(meas.K,1);
est.L= cell(meas.K,1);

%filter parameters
filter.T_max= 100;                  %maximum number of tracks
filter.track_threshold= 1e-3;       %threshold to prune tracks

filter.H_bth= 5;                    %requested number of birth components/hypotheses (for LMB to GLMB casting before update)
filter.H_sur= 200;                  %requested number of surviving components/hypotheses (for LMB to GLMB casting before update)
filter.H_upd= 200;                  %requested number of updated components/hypotheses (for GLMB update)
filter.H_max= 200;                  %cap on number of posterior components/hypotheses (not used yet)
filter.hyp_threshold= 1e-12;         %pruning threshold for components/hypotheses (not used yet)

filter.npt= 1000;                   %number of particles per track
filter.nth= 500;                    %threshold on effective number of particles before resampling (not used here, resampling is forced at every step, otherwise number of particles per track grows)

filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output

est.filter= filter;

%=== Filtering

%initial prior
tt_lmb_update= cell(0,1);      %track table for LMB (cell array of structs for individual tracks)

%recursive filtering
for k=1:meas.K
    
    %prediction
    [tt_lmb_birth,tt_lmb_survive]= lmbpredict(tt_lmb_update,model,filter,k);            T_predict= length(tt_lmb_birth)+length(tt_lmb_survive);

    %update
    glmb_predict= castlmbpred(tt_lmb_birth,tt_lmb_survive,filter);
    glmb_update= update(glmb_predict,model,filter,meas,k);
    tt_lmb_update= glmb2lmb(glmb_update);                                               T_posterior= length(tt_lmb_update);
       
    %pruning, truncation and track cleanup
    tt_lmb_update= clean_lmb(tt_lmb_update,filter);                                     T_clean= length(tt_lmb_update);
    
    %state estimation
    [est.X{k},est.N(k),est.L{k}]= extract_estimates(tt_lmb_update,model);
    
    %display diagnostics
    display_diaginfo(tt_lmb_update,k,est,filter,T_predict,T_posterior,T_clean);
    
end
end

function [tt_lmb_birth,tt_lmb_survive]= lmbpredict(tt_lmb_update,model,filter,k)
%---generate birth tracks
tt_lmb_birth= cell(length(model.r_birth),1);                                                                            %initialize cell array
for tabbidx=1:length(model.r_birth)
    tt_lmb_birth{tabbidx}.r= model.r_birth(tabbidx);                                                                    %birth prob for birth track
    tt_lmb_birth{tabbidx}.x= gen_gms(model.w_birth{tabbidx},model.m_birth{tabbidx},model.P_birth{tabbidx},filter.npt);  %samples for birth track
    tt_lmb_birth{tabbidx}.w= ones(filter.npt,1)/filter.npt;                                                             %weights of samples for birth track
    tt_lmb_birth{tabbidx}.l= [k;tabbidx];                                                                               %track label
end

%---generate surviving tracks
tt_lmb_survive= cell(length(tt_lmb_update),1);                                                                              %initialize cell array
for tabsidx=1:length(tt_lmb_update)
    wtemp_predict= compute_pS(model,tt_lmb_update{tabsidx}.x).*tt_lmb_update{tabsidx}.w(:); xtemp_predict= gen_newstate_fn(model,tt_lmb_update{tabsidx}.x,'noise');      %particle prediction
    tt_lmb_survive{tabsidx}.r= sum(wtemp_predict)*tt_lmb_update{tabsidx}.r;                                                 %predicted existence probability for surviving track
    tt_lmb_survive{tabsidx}.x= xtemp_predict;                                                                               %samples for surviving track
    tt_lmb_survive{tabsidx}.w= wtemp_predict/sum(wtemp_predict);                                                            %weights of samples for predicted track
    tt_lmb_survive{tabsidx}.l= tt_lmb_update{tabsidx}.l;                                                                    %track label
end
end



function glmb_predict= castlmbpred(tt_lmb_birth,tt_lmb_survive,filter)
%express birth and surviving LMBs in GLMB structure
glmb_birth= lmb2glmb(tt_lmb_birth,filter.H_bth);
glmb_survive= lmb2glmb(tt_lmb_survive,filter.H_sur);

%generate predicted hypotheses/components (by convolution of birth and survive GLMBs)
glmb_predict.tt= cat(1,glmb_birth.tt,glmb_survive.tt);
for bidx= 1:length(glmb_birth.w)
    for sidx= 1:length(glmb_survive.w)
        hidx= (bidx-1)*length(glmb_survive.w)+sidx;
        glmb_predict.w(hidx)= glmb_birth.w(bidx)*glmb_survive.w(sidx);
        glmb_predict.I{hidx}= [glmb_birth.I{bidx}; length(glmb_birth.tt)+glmb_survive.I{sidx}];
        glmb_predict.n(hidx)= glmb_birth.n(bidx)+glmb_survive.n(sidx);
    end
end
glmb_predict.w= glmb_predict.w/sum(glmb_predict.w);

%extract cardinality distribution
for card=0:max(glmb_predict.n)
    glmb_predict.cdn(card+1)= sum(glmb_predict.w(glmb_predict.n==card));
end
end



function glmb= lmb2glmb(tt_lmb,H_req)
%express LMB in GLMB structure
rvect= get_rvals(tt_lmb);                                   %vector of existence probabilities
costv= rvect./(1-rvect);                                    %cost vector
neglogcostv= -log(costv);                                   %negative log cost
[paths,nlcost]= kshortestwrap_pred(neglogcostv,H_req);      %k-shortest path to calculate k-best surviving hypotheses/components
glmb.tt= tt_lmb;
for hidx=1:length(nlcost)
    glmb.w(hidx)= sum(log(1-rvect))-nlcost(hidx);           %weight of hypothesis/component
    glmb.I{hidx}= paths{hidx}(:);                           %tracks in hypothesis/component
    glmb.n(hidx)= length(paths{hidx}(:));                   %cardinality of hypothesis/component
end
glmb.w= exp(glmb.w-logsumexp(glmb.w));                      %normalize weights

%extract cardinality distribution
for card=0:max(glmb.n)
    glmb.cdn(card+1)= sum(glmb.w(glmb.n==card));            %extract probability of n targets
end
end



function glmb_update= update(glmb_predict,model,filter,meas,k)
%create updated tracks (single target Bayes update)
m= size(meas.Z{k},2);                                   %number of measurements
tt_update= cell((1+m)*length(glmb_predict.tt),1);       %initialize cell array
%missed detection tracks (legacy tracks)
for tabidx= 1:length(glmb_predict.tt)
    tt_update{tabidx}= glmb_predict.tt{tabidx};         %same track table
end

%measurement updated tracks (all pairs)
allcostm= zeros(length(glmb_predict.tt),m);                                                 %global cost matrix 
for emm= 1:m
    for tabidx= 1:length(glmb_predict.tt)
        stoidx= length(glmb_predict.tt)*emm + tabidx; %index of predicted track i updated with measurement j is (number_predicted_tracks*j + i)
        w_temp= compute_pD(model,glmb_predict.tt{tabidx}.x).*glmb_predict.tt{tabidx}.w(:).*compute_likelihood(model,meas.Z{k}(:,emm),glmb_predict.tt{tabidx}.x)'; x_temp= glmb_predict.tt{tabidx}.x;  %weight update for this track and this measuremnent      
        tt_update{stoidx}.x= x_temp;                                                        %particles for updated track
        tt_update{stoidx}.w= w_temp/sum(w_temp);                                            %weights of partcles for updated track
        tt_update{stoidx}.l = glmb_predict.tt{tabidx}.l;                                    %track label
        allcostm(tabidx,emm)= sum(w_temp);                                                  %predictive likelihood
    end
end
glmb_update.tt= tt_update;                                                                  %copy track table back to GLMB struct

%precalculation loop for average detection/missed probabilities
avpd= zeros(length(glmb_predict.tt),1);
for tabidx=1:length(glmb_predict.tt)
    avpd(tabidx)= glmb_predict.tt{tabidx}.w(:)'*compute_pD(model,glmb_predict.tt{tabidx}.x)+eps(0);
end
avqd= 1-avpd;

%component updates
if m==0 %no measurements means all missed detections
    for pidx=1:length(glmb_predict.w)
        glmb_update.w(pidx)= -model.lambda_c+sum(log(avqd(glmb_predict.I{pidx})))+log(glmb_predict.w(pidx));            %hypothesis/component weight
    end
    glmb_update.I= glmb_predict.I;                                                                                      %hypothesis/component tracks (via indices to track table)
    glmb_update.n= glmb_predict.n;                                                                                      %hypothesis/component cardinality
else %loop over predicted components/hypotheses
    runidx= 1;
    for pidx=1:length(glmb_predict.w)
        if glmb_predict.n(pidx)==0 %no target means all clutter
            glmb_update.w(runidx)= -model.lambda_c+m*log(model.lambda_c*model.pdf_c)+log(glmb_predict.w(pidx));
            glmb_update.I{runidx}= glmb_predict.I{pidx};
            glmb_update.n(runidx)= glmb_predict.n(pidx);
            runidx= runidx+1;
        else %otherwise perform update for component
            %calculate best updated hypotheses/components
            costm= allcostm(glmb_predict.I{pidx},:)./(model.lambda_c*model.pdf_c*repmat(avqd(glmb_predict.I{pidx}),[1 m]));                     %cost matrix
            neglogcostm= -log(costm);                                                                                                           %negative log cost
            [uasses,nlcost]= mbestwrap_updt_custom(neglogcostm,round(filter.H_upd*sqrt(glmb_predict.w(pidx))/sum(sqrt(glmb_predict.w))));       %murty's algo to calculate m-best assignment hypotheses/components
            
            %generate corrresponding surviving hypotheses/components
            for hidx=1:length(nlcost)
                update_hypcmp_tmp= uasses(hidx,:)';
                glmb_update.w(runidx)= -model.lambda_c+m*log(model.lambda_c*model.pdf_c)+sum(log(avqd(glmb_predict.I{pidx})))+log(glmb_predict.w(pidx))-nlcost(hidx);       %hypothesis/component weight
                glmb_update.I{runidx}= length(glmb_predict.tt).*update_hypcmp_tmp+glmb_predict.I{pidx};                                                                     %hypothesis/component tracks (via indices to track table)
                glmb_update.n(runidx)= glmb_predict.n(pidx);                                                                                                                %hypothesis/component cardinality
                runidx= runidx+1;
            end
        end
    end
end
glmb_update.w= exp(glmb_update.w-logsumexp(glmb_update.w));                                                                                                                 %normalize weights

%extract cardinality distribution
for card=0:max(glmb_update.n)
    glmb_update.cdn(card+1)= sum(glmb_update.w(glmb_update.n==card));                                                                                                       %extract probability of n targets
end
end



function tt_lmb= glmb2lmb(glmb)

%find unique labels (with different possibly different association histories)
lmat= zeros(2,length(glmb.tt),1);
for tabidx= 1:length(glmb.tt)
    lmat(:,tabidx)= glmb.tt{tabidx}.l;
end
lmat= lmat';

[cu,~,ic]= unique(lmat,'rows'); cu= cu';

%initialize LMB struct
tt_lmb= cell(size(cu,2),1);
for tabidx=1:length(tt_lmb)
   tt_lmb{tabidx}.r= 0;
   tt_lmb{tabidx}.x= [];
   tt_lmb{tabidx}.w= [];
   tt_lmb{tabidx}.l= cu(:,tabidx);
end

%extract individual tracks
for hidx=1:length(glmb.w)
   for t= 1:glmb.n(hidx)
      trkidx= glmb.I{hidx}(t);
      newidx= ic(trkidx);
      tt_lmb{newidx}.x= cat(2,tt_lmb{newidx}.x,glmb.tt{trkidx}.x);
      tt_lmb{newidx}.w= cat(1,tt_lmb{newidx}.w,glmb.w(hidx)*glmb.tt{trkidx}.w);
   end
end

%extract existence probabilities and normalize track weights
for tabidx=1:length(tt_lmb)
   tt_lmb{tabidx}.r= sum(tt_lmb{tabidx}.w);
   tt_lmb{tabidx}.w= tt_lmb{tabidx}.w/tt_lmb{tabidx}.r;
end

end



function tt_lmb_out= clean_lmb(tt_lmb_in,filter)
%prune tracks with low existence probabilities
rvect= get_rvals(tt_lmb_in);
idxkeep= find(rvect > filter.track_threshold);
rvect= rvect(idxkeep);
tt_lmb_out= tt_lmb_in(idxkeep);

%enforce cap on maximum number of tracks
if length(tt_lmb_out) > filter.T_max
    [~,idxkeep]= sort(rvect,'descend');
    tt_lmb_out= tt_lmb_out(idxkeep);   
end

%cleanup tracks
for tabidx=1:length(tt_lmb_out)
    xtemptemp= tt_lmb_out{tabidx}.x;
    wtemptemp= tt_lmb_out{tabidx}.w;
    rspidx= randsample(length(wtemptemp),filter.npt,true,wtemptemp); %rspidx= resample(wtemptemp,filter.npt);
    tt_lmb_out{tabidx}.x= xtemptemp(:,rspidx);
    tt_lmb_out{tabidx}.w= ones(filter.npt,1)'/filter.npt;
end
end



function rvect= get_rvals(tt_lmb)                           %function to extract vector of existence probabilities from LMB track table
rvect= zeros(length(tt_lmb),1);
for tabidx=1:length(tt_lmb)
   rvect(tabidx)= tt_lmb{tabidx}.r; 
end
end



function [X,N,L]=extract_estimates(tt_lmb,model)
%extract estimates via MAP cardinality and corresponding tracks
rvect= get_rvals(tt_lmb);
cdn= prod(1-rvect)*esf(rvect./(1-rvect));
[~,mode] = max(cdn);
N = min(length(rvect),mode-1);
X= zeros(model.x_dim,N);
L= zeros(2,N);

[~,idxcmp]= sort(rvect,'descend');
for n=1:N
    [~,idxtrk]= max(tt_lmb{idxcmp(n)}.w);
    X(:,n)= tt_lmb{idxcmp(n)}.x*tt_lmb{idxcmp(n)}.w(:);
    L(:,n)= tt_lmb{idxcmp(n)}.l;
end
end



function display_diaginfo(tt_lmb,k,est,filter,T_predict,T_posterior,T_clean)
rvect= get_rvals(tt_lmb);
cdn= prod(1-rvect)*esf(rvect./(1-rvect));
eap= (0:(length(cdn)-1))*cdn(:);
var= (0:(length(cdn)-1)).^2*cdn(:)-((0:(length(cdn)-1))*cdn(:))^2;
if ~strcmp(filter.run_flag,'silence')
    disp([' time= ',num2str(k),...
        ' #eap cdn=' num2str(eap),...
        ' #var cdn=' num2str(var,4),...
        ' #est card=' num2str(est.N(k),4),...
        ' #trax pred=' num2str(T_predict,4),...
        ' #trax post=' num2str(T_posterior,4),...
        ' #trax updt=',num2str(T_clean,4)   ]);
end
end