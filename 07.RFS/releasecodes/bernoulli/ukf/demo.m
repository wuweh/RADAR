% This is a demo script for the Bernoulli filter with RFS observations proposed in
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

model= gen_model;
truth= gen_truth(model);
meas=  gen_meas(model,truth);
est=   run_filter(model,meas);
handles= plot_results(model,truth,meas,est);