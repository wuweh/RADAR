% This is a demo script for the SMC-PHD filter proposed in
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

model= gen_model;
truth= gen_truth(model);
meas=  gen_meas(model,truth);
est=   run_filter(model,meas);
handles= plot_results(model,truth,meas,est);