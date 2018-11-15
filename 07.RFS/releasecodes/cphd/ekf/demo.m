% This is a demo script for the GMCPHD filter proposed in
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

model= gen_model;
truth= gen_truth(model);
meas=  gen_meas(model,truth);
est=   run_filter(model,meas);
handles= plot_results(model,truth,meas,est);