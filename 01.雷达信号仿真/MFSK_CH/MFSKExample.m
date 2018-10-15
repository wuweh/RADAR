%目前存在的问题：
%           1) 由于测量误差以及系统的误差，导致接收信号的相位差不精确（这应该是最致命的问题），从而影响目标信息的结算
%           2) 会出现不同的距离和速度的目标，在频率谱上无法分辨的现象（暂时不明白）
%           3）实际应用中，需要独立采样MFSK的两组波形，对于软件流程控制有一定的要求

clc;
close all;
clear all;

%生成波形参数
[fmcwwaveform,target,tgtmotion,channel,transmitter,receiver,...
    sensormotion,c,fc,lambda,fs,maxbeatfreq] = helperMFSKSystemSetup;

%设置波形参数
mfskwaveform =  phased.MFSKWaveform(...
                'SampleRate',151e6,...
                'SweepBandwidth',150e6,...
                'StepTime',2e-6,...
                'StepsPerSweep',1024,...
                'FrequencyOffset',-294e3,...
                'OutputFormat','Sweeps',...
                'NumSweeps',1);

%%
% The figure below shows the spectrogram of the waveform. It is zoomed into
% a small interval to better reveal the time frequency characteristics of
% the waveform.

%显示波形频谱
numsamp_step = round(mfskwaveform.SampleRate*mfskwaveform.StepTime);
sig_display = mfskwaveform();
% spectrogram(sig_display(1:8192),kaiser(3*numsamp_step,100),...
%     ceil(2*numsamp_step),linspace(0,4e6,2048),mfskwaveform.SampleRate,...
%     'yaxis','reassigned','minthreshold',-60)

%%
% Next, simulate the return of the system. Again, only 1 sweep is needed to
% estimate the range and Doppler.

Nsweep = 1;
release(channel);
channel.SampleRate = mfskwaveform.SampleRate;
release(receiver);
receiver.SampleRate = mfskwaveform.SampleRate;

xr = helperFMCWSimulate(Nsweep,mfskwaveform,sensormotion,tgtmotion,...
    transmitter,channel,target,receiver);

%%
% The subsequent processing samples the return echo at the end of each step
% and group the sampled signals into two sequences corresponding to two
% sweeps. Note that the sampling frequency of the resulting sequence is now
% proportional to the time at each step, which is much less compared the
% original sample rate.

x_dechirp = reshape(xr(numsamp_step:numsamp_step:end),2,[]).';
fs_dechirp = 1/(2*mfskwaveform.StepTime);


%%
% As in the case of FMCW signals, the MFSK waveform is processed in the
% frequency domain. Next figures shows the frequency spectrums of the
% received echos corresponding to the two sweeps.

xf_dechirp = fft(x_dechirp);
num_xf_samp = size(xf_dechirp,1);
beatfreq_vec = (0:num_xf_samp-1).'/num_xf_samp*fs_dechirp;

figure;
plot(beatfreq_vec(1:100)/1e3,abs(xf_dechirp(1:100,1)),'-b*');hold on;grid on;
plot(beatfreq_vec(1:100)/1e3,abs(xf_dechirp(1:100,2)),'-r*');hold on;grid on;
ylabel('Magnitude');
title('Frequency spectrum for sweep');
xlabel('Frequency (kHz)')

%%
% Note that there are two peaks in each frequency spectrum indicating two
% targets. In addition, the peaks are at the identical locations in both
% returns so there is no ghost targets.
%
% To detect the peaks, one can use a CFAR detector. Once detected, the beat
% frequencies as well as the phase differences between two spectra are
% computed at the peak locations.

% cfar = phased.CFARDetector('ProbabilityFalseAlarm',1e-1,...
%     'NumTrainingCells',8);
% 
% peakidx = cfar(abs(xf_dechirp(:,1)),1:num_xf_samp);
% Fbeat = beatfreq_vec(peakidx);
% phi = angle(xf_dechirp(peakidx_1,2))-angle(xf_dechirp(peakidx_1,1));

%目标检测
k = 1;
for n= 10:1:512/2
    if (abs(xf_dechirp(n,1))>abs(xf_dechirp(n-1,1))) && ...
       (abs(xf_dechirp(n,1))>abs(xf_dechirp(n-2,1))) && ...
       (abs(xf_dechirp(n,1))>abs(xf_dechirp(n-3,1))) && ...
       (abs(xf_dechirp(n,1))>abs(xf_dechirp(n-4,1))) && ...
       (abs(xf_dechirp(n,1))>abs(xf_dechirp(n-5,1))) && ...
       (abs(xf_dechirp(n,1))>abs(xf_dechirp(n+1,1))) && ...
       (abs(xf_dechirp(n,1))>abs(xf_dechirp(n+2,1))) && ...
       (abs(xf_dechirp(n,1))>abs(xf_dechirp(n+3,1))) && ...
       (abs(xf_dechirp(n,1))>abs(xf_dechirp(n+4,1))) && ...
       (abs(xf_dechirp(n,1))>abs(xf_dechirp(n+5,1))) && ...
       (abs(xf_dechirp(n,1))>0.02)
        Fbeat(k) = beatfreq_vec(n);
        peakidx_1(k) = n;
        k = k + 1;
    end
end
Fbeat = Fbeat'

%计算相位差
phi = angle(xf_dechirp(peakidx_1,2))-angle(xf_dechirp(peakidx_1,1));


%%
% Finally, the beat frequencies and phase differences are used to estimate
% the range and speed. Depending on how one constructs the phase
% difference, the equations are slightly different. For the approach shown
% in this example, it can be shown that the range and speed satisfies the
% following relation:
%
% $$f_b = -\frac{2v}{\lambda}+ \frac{2\beta R}{c} $$
%
% $$\Delta\phi = -\frac{4\pi T_sv}{\lambda} + \frac{4\pi f_{offset}R}{c}$$
%
% where $f_b$ is the beat frequency, $\Delta\phi$ is the phase difference,
% $\lambda$ is the wavelength, $c$ is the propagation speed, $T_s$ is the
% step time, $f_{offset}$ is the frequency offset, $\beta$ is the sweep
% slope, $R$ is the range, and $v$ is the speed. Based on the equation, the
% range and speed are estimated below:

%解算目标信息（距离、速度）
sweep_slope = mfskwaveform.SweepBandwidth/...
    (mfskwaveform.StepsPerSweep*mfskwaveform.StepTime);
temp =  [1 sweep_slope;mfskwaveform.StepTime mfskwaveform.FrequencyOffset]\...
        [Fbeat phi/(2*pi)].';

r_est = c*temp(2,:)/2
v_est = lambda*temp(1,:)/(-2)*3.6



