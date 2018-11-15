% This is a demo script for the CBMeMBer filter proposed in
% (without track labelling)
% B.-T. Vo, B.-N. Vo, and A. Cantoni, "The Cardinality Balanced Multi-target Multi-Bernoulli filter and its implementations," IEEE Trans. Signal Processing, Vol. 57, No. 2, pp. 409–423, 2009. 
% http://ba-ngu.vo-au.com/vo/VVCmemberSP09.pdf
% ---BibTeX entry
% @ARTICLE{CBMEMBER,
% author={B.-T. Vo and B.-N. Vo and A. Cantoni},
% journal={IEEE Transactions on Signal Processing},
% title={The Cardinality Balanced Multi-Target Multi-Bernoulli Filter and Its Implementations},
% year={2009},
% month={Feb},
% volume={57},
% number={2},
% pages={409-423}} 
%---

model= gen_model;
truth= gen_truth(model);
meas=  gen_meas(model,truth);
est=   run_filter(model,meas);
handles= plot_results(model,truth,meas,est);