% This is a demo script for the SMC implementation of the CPHD filter given in
% (assuming Poisson clutter)
% B.-T. Vo, B.-N. Vo and A. Cantoni, "Analytic implementations of the Cardinalized Probability Hypothesis Density Filter," IEEE Trans Signal Processing, Vol. 55,  No. 7, part 2,  pp. 3553-3567, 2007.
% http://ba-ngu.vo-au.com/vo/VVC_CPHD_SP07.pdf
% ---BibTeX entry
% @ARTICLE{GMCPHD,
% author={B.-T. Vo and B.-N. Vo and A. Cantoni},
% journal={IEEE Transactions on Signal Processing},
% title={Analytic Implementations of the Cardinalized Probability Hypothesis Density Filter},
% year={2007},
% month={July},
% volume={55},
% number={7},
% pages={3553-3567}} 
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