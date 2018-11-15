% This is a demo script for the Generalized Labeled Multi-Bernoulli filter proposed in
% B.-T. Vo, and B.-N. Vo, "Labeled Random Finite Sets and Multi-Object Conjugate Priors," IEEE Trans. Signal Processing, Vol. 61, No. 13, pp. 3460-3475, 2013.
% http://ba-ngu.vo-au.com/vo/VV_Conjugate_TSP13.pdf
% with corresponding implementation details given in
% B.-N. Vo, B.-T. Vo, and D. Phung, "Labeled Random Finite Sets and the Bayes Multi-Target Tracking Filter," IEEE Trans. Signal Processing, Vol. 62, No. 24, pp. 6554-6567, 2014
% http://ba-ngu.vo-au.com/vo/VVP_GLMB_TSP14.pdf
%
% Note 1: no lookahead PHD/CPHD allocation is implemented in this code, a simple proportional weighting scheme is used for readability
% Note 2: the simple example used here is the same as in the CB-MeMBer filter code for a quick demonstration and comparison purposes
% Note 3: more difficult scenarios require a better lookahead and more components/hypotheses (and exec time)
% ---BibTeX entry
% @ARTICLE{GLMB1,
% author={B.-T. Vo and B.-N. Vo
% journal={IEEE Transactions on Signal Processing},
% title={Labeled Random Finite Sets and Multi-Object Conjugate Priors},
% year={2013},
% month={Jul}
% volume={61},
% number={13},
% pages={3460-3475}}
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

model= gen_model;
truth= gen_truth(model);
meas=  gen_meas(model,truth);
est=   run_filter(model,meas);
handles= plot_results(model,truth,meas,est);