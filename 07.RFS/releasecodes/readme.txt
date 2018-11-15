====================[http://ba-tuong.vo-au.com/codes.html]====================

This is a beta release for a suite of MATLAB based RFS filtering/tracking codes.

The "_common" subdirectory of shared functions should be added to your MATLAB path.

The remaining directories (except for "robust") contain various filters:
	"singletarget"		Single target filter with RFS observations 
	"bernoulli" 		Bernoulli filter
	"phd"			PHD filter
	"cphd"			CPHD filter
	"cbmember"		CBMeMBer filter
	"glmb"			Generalized Labeled Multi-Bernoulli Filter
	"lmb"			Labeled Multi-Bernoulli filter

Each filter has different implementation folders:
	"gms"			Gaussian mixture solution for linear Gaussian models
	"ekf"			EKF approximation for non-linear models
	"ukf"			UKF approximation for non-linear models
	"smc"			SMC implementation for nonlinear models

The "robust" directory contains implementations for
	"lccphd"		Unknown lambda-CPHD filter (with GM/EKF/UKF/SMC versions)
	"pdcphd"		Unknown P_D CPHD filter	(with BGM/EKF/UKF/SMC versions)
	"jointcphd"		Jointly unknown lambda and P_D CPHD filter (with BGM/EKF/UKF/SMC versions)
	"jointcbmember"		Jointly unknown lambda and P_D CBMeMBer filter (SMC version only; but unlike the CPHD variants which estimate the clutter rate only (with known map), this estimates clutter maps on the fly)

For each filter and within each implementation folder, run the "demo" script to see a preconfigured example.
e.g. to run the CPHD filter with a UKF based implementation:

>> addpath _common
>> cd cphd/ukf/
>> demo

Notes:

1. All linear Gaussian examples use a 4D CV model (x/y position and velocity) with 2D observations (position only).

2. All non linear examples use a 5D CT model (x/y position, velocity and unknown turn rate) with 2D observations (range and bearing).

3. The demo scripts all use the same basic calling sequence: 

model= gen_model;					%generate model parameters
truth= gen_truth(model);				%generate ground truths
meas=  gen_meas(model,truth);				%generate measurements
est=   run_filter(model,meas);				%perform filtering
handles= plot_results(model,truth,meas,est);		%plotting

4. These codes are provided for academic/research purposes only. They are designed for readability and demonstration thus they are not optimized for speed. The authors do not provide any guarantees implicit or otherwise regarding these codes. For commercial grade codes please contact the authors directly (ba-tuong.vo@curtin.edu.au, ba-ngu.vo@curtin.edu.au). The authors acknowledge the input of and discussions with Dr Michael Beard and Dr Stephan Reuter respectively in the development of these particular GLMB and LMB codes. The codes for "robust" filtering were developed by Dr Du Yong Kim (duyongkim@gmail.com/duyong.kim@curtin.edu.au).

