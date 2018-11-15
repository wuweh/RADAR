% This a demo script for the SMC implementation of the lambda-CPHD filter proposed in
% R. Mahler, B.-T. Vo and B.-N. Vo, "CPHD Filtering with unknown clutter rate and detection profile," IEEE Transactions on Signal Processing, Vol. 59, No. 8, pp. 3497-3513, 2011.
% http://ba-ngu.vo-au.com/vo/MVV_PHDrobust.pdf

% There are three versions of filters for this paper
% 1) Lambda-CPHD, 2) pD-CPHD, and 3) Lambda-pD-CPHD as sequentially described in the paper.
% This is the code for Lambda-CPHD

% ---BibTeX entry
% @ARTICLE{RobustCPHD,
% author={R. P. S. Mahler and B.-T. Vo and B.-N. Vo},
% journal={IEEE Transactions on Signal Processing},
% title={CPHD Filtering With Unknown Clutter Rate and Detection Profile},
% year={2011},
% month={August},
% volume={59},
% number={8},
% pages={3497-3513}} 
%---
% 
% based on the SMC implementation of the PHD filter given in 
%
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


model= gen_model;
truth= gen_truth(model);
meas=  gen_meas(model,truth);
est=   run_filter(model,meas);
handles= plot_results(model,truth,meas,est);