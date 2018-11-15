% This is a demo script for the single target filter with RFS observations proposed in
% (assuming Poisson clutter and no extraneous measurements)
% B.-T. Vo, B.-N. Vo, and A. Cantoni, "Bayesian filtering with random finite set observations," IEEE Trans. Signal Processing, Vol. 56, No. 4, pp. 1313-1326, 2008.
% http://ba-ngu.vo-au.com/vo/VVCsingletargetSP08.pdf
% ---BibTeX entry
% @ARTICLE{STF, 
% author={B.-T.Vo and B.-N. Vo and A. Cantoni},
% journal={IEEE Transactions on Signal Processing},
% title={Bayesian Filtering With Random Finite Set Observations},
% year={2008},
% month={April},
% volume={56},
% number={4},
% pages={1313-1326}}
%---
clc;clear all; close all;
model= gen_model;               %CTģ�͡����û���������Q��R�ȵȣ�
truth= gen_truth(model);        %����Ŀ���˶��켣
meas=  gen_meas(model,truth);   %��������㼣

est=   run_filter(model,meas);          

handles= plot_results(model,truth,meas,est);  %��ͼ