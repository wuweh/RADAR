% This is a demo script for the GM-PHD filter proposed in
% (assuming no target spawning)
% B.-N. Vo, and W. K. Ma, "The Gaussian mixture Probability Hypothesis Density Filter," IEEE Trans Signal Processing, Vol. 54, No. 11, pp. 4091-4104, 2006.
% http://ba-ngu.vo-au.com/vo/VM_GMPHD_SP06.pdf
% ---BibTeX entry
% @ARTICLE{GMPHD,
% author={B.-N. Vo and W.-K. Ma},
% journal={IEEE Transactions on Signal Processing},
% title={The Gaussian Mixture Probability Hypothesis Density Filter},
% year={2006},
% month={Nov},
% volume={54},
% number={11},
% pages={4091-4104}}
%---

model= gen_model;
truth= gen_truth(model);
meas=  gen_meas(model,truth);
est=   run_filter(model,meas);
handles= plot_results(model,truth,meas,est);