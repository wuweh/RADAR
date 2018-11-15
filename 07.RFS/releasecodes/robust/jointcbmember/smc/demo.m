% This is the demo script for the Robust CBMeMBer filter proposed in
% (without track labelling)
% B.-T. Vo, B.-N. Vo, R. Hoseinnezhad, and R. Mahler "Robust Multi-Bernoulli Filtering," IEEE Journal on Selected Topics in Signal Processing, Vol. 7, No. 3, pp. 399-409, 2013.
% http://ba-ngu.vo-au.com/vo/VVHM_JSSP13.pdf
% ---BibTeX entry
% @ARTICLE{RobustCBMEMBER,
% author={B.-T. Vo and B.-N. Vo and R. Hoseinnezhad and R. Mahler},
% journal={IEEE Transactions on Signal Processing},
% title={Robust Multi-Bernoulli Filtering},
% year={2013},
% month={Jun},
% volume={7},
% number={3},
% pages={399-409}} 
%---

model= gen_model;
truth= gen_truth(model);
meas=  gen_meas(model,truth);
est=   run_filter(model,meas);
handles= plot_results(model,truth,meas,est);