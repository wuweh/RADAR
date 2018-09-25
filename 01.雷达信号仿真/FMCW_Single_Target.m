clc;
clear all;
close all;

fc = 77e9;  %77GHz
c = 3e8;      %光速
lambda = c/fc;  %波长

%%
% The sweep time can be computed based on the time needed for the signal to
% travel the unambiguous maximum range. In general, for an FMCW radar
% system, the sweep time should be at least 5 to 6 times the round trip
% time. 
% This example uses a factor of 5.5.
% The function computes 2*r/c.
%由最远检测距离确定扫频时间
range_max = 200;                                        %最远检测距离
tm = 5.5*range2time(range_max,c);           %扫频时间

%%
% The sweep bandwidth can be determined according to the range resolution
% and the sweep slope is calculated using both sweep bandwidth and sweep
% time.

% 距离分辨率（range resolution) = c * Tp / 2 
% range resolution = C/(2*B);
% 带宽 B 近似等于 1/Tp, 所以带宽越大，分辨率越高。
range_res = 1;
bw = range2bw(range_res,c);    %由距离分辨率确定扫频带宽
sweep_slope = bw/tm;                %扫频斜率 = 扫频带宽/扫频时间

%%
% Because an FMCW signal often occupies a huge bandwidth, setting the
% sample rate blindly to twice the bandwidth often stresses the capability
% of A/D converter hardware. To address this issue, one can often choose a
% lower sample rate. Two things can be considered here:
%
% # For a complex sampled signal, the sample rate can be set to the same as
% the bandwidth.
% # FMCW radars estimate the target range using the beat frequency embedded
% in the dechirped signal. The maximum beat frequency the radar needs to
% detect is the sum of the beat frequency corresponding to the maximum
% range and the maximum Doppler frequency. Hence, the sample rate only
% needs to be twice the maximum beat frequency.
%
% In this example, the beat frequency corresponding to the maximum range is
% given by

% The function computes 2*range_max*sweep_slope/c 
%由最远检测距离带来的最大差频频率
fr_max = range2beat(range_max,sweep_slope,c);

%%
% In addition, the maximum speed of a traveling car is about 230 km/h.
% Hence the maximum Doppler shift and the maximum beat frequency can be
% computed as
%由最大速度引起的多普勒频率：fd_max = 2*v_max/lambda;
%2*v_max是因为两车可能是相背行驶
v_max = 230*1000/3600;
fd_max = speed2dop(2*v_max,lambda);

%由此可以得出最终的最大差频频率
fb_max = fr_max+fd_max;


%%
% This example adopts a sample rate of the larger of twice the maximum beat
% frequency and the bandwidth.

%确定采样频率
fs = max(2*fb_max,bw);

%%
% The following table summarizes the radar parameters.
% 
%  System parameters            Value
%  ----------------------------------
%  Operating frequency (GHz)    77
%  Maximum target range (m)     200
%  Range resolution (m)         1
%  Maximum target speed (km/h)  230
%  Sweep time (microseconds)    7.33
%  Sweep bandwidth (MHz)        150
%  Maximum beat frequency (MHz) 27.30
%  Sample rate (MHz)            150


%%
% With all the information above, one can set up the FMCW waveform used
% in the radar system.

%生成FMCW波形参数
hwav = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw, 'SampleRate',fs);

%%
% This is a up-sweep linear FMCW signal, often referred to as sawtooth
% shape. One can examine the time-frequency plot of the generated signal.

% 形成波形
s = step(hwav);
%实部即信号的幅度！！！！
subplot(211); plot(0:1/fs:tm-1/fs,real(s));
xlabel('Time (s)'); ylabel('Amplitude (v)');
title('FMCW signal'); axis tight;

%spectrogram: 用短时傅里叶变换
%16阶hanmming windows；
%0 sample overlap
% 4096点FFT
%采样频率为fs
%发射波形频谱图
subplot(212); spectrogram(s,32,16,4096,fs,'yaxis');
title('FMCW signal spectrogram');

%% Target Model
% The target of an ACC radar is usually a car in front of it. This example
% assumes the target car is moving 50 m ahead of the car with the
% radar, at a speed of 96 km/h along the x-axis. 
%
% The radar cross section of a car, according to [1], can be computed based
% on the distance between the radar and the target car.

car_dist = 30;
car_speed = 170*1000/3600;

%目标车辆的RCS值
car_rcs = db2pow(min(10*log10(car_dist)+5,20));

%c指的是光速，fc指的是载波频率
hcar = phased.RadarTarget('MeanRCS',car_rcs,'PropagationSpeed',c,...
    'OperatingFrequency',fc);

%定义目标车辆的初始位置[x,y,z]和初始速度[vx,vy,vz]
hcarplatform = phased.Platform('InitialPosition',[car_dist;5;0.5],...
    'Velocity',[car_speed;0;0]);

%%
% The propagation model is assumed to be free space.

hchannel = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

%% Radar System Setup
% The rest of the radar system includes the transmitter, the receiver, and
% the antenna. This example uses the parameters presented in [1]. Note that
% this example models only main components and omits the effect from other
% components, such as coupler and mixer. In addition, for the sake of
% simplicity, the antenna is assumed to be isotropic and the gain of the
% antenna is included in the transmitter and the receiver.

%t天线孔径大小 m^2
ant_aperture = 6.06e-4;                         % in square meter
%转换为天线增益
ant_gain = aperture2gain(ant_aperture,lambda);  % in dB

tx_ppower = db2pow(5)*1e-3;                     % in watts
tx_gain = 9+ant_gain;                           % in dB

rx_gain = 15+ant_gain;                          % in dB
rx_nf = 4.5;                                    % in dB

htx = phased.Transmitter('PeakPower',tx_ppower,'Gain',tx_gain);
hrx = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,...
    'SampleRate',fs);

%%
% Automotive radars are generally mounted on a vehicles, so they are often
% in motion. This example assumes the radar is traveling at a speed of 100
% km/h along x-axis. So the target car is approaching the radar at a
% relative speed of 4 km/h.

%本车速度 100km/h和初始位置信息
radar_speed = 100*1000/3600;
hradarplatform = phased.Platform('InitialPosition',[0;0;0.5],...
    'Velocity',[radar_speed;0;0]);

%% Radar Signal Simulation
% As briefly mentioned in earlier sections, an FMCW radar measures the
% range by examining the beat frequency in the dechirped signal. To extract
% this frequency, a dechirp operation is performed by mixing the received
% signal with the transmitted signal. After the mixing, the dechirped
% signal contains only individual frequency components that correspond to
% the target range.
% 
% In addition, even though it is possible to extract the Doppler
% information from a single sweep, the Doppler shift is often extracted
% among several sweeps because within one pulse, the Doppler frequency is
% indistinguishable from the beat frequency. To measure the range and
% Doppler, an FMCW radar typically performs the following operations:
%
% # The waveform generator generates the FMCW signal.
% # The transmitter and the antenna amplify the signal and radiate the
% signal into space.
% # The signal propagates to the target, gets reflected by the target, and
% travels back to the radar.
% # The receiving antenna collects the signal.
% # The received signal is dechirped and saved in a buffer.
% # Once a certain number of sweeps fill the buffer, the Fourier transform
% is performed in both range and Doppler to extract the beat frequency as
% well as the Doppler shift. One can then estimate the range and speed of
% the target using these results. Range and Doppler can also be shown as an
% image and give an intuitive indication of where the target is in the
% range and speed domain.
%
% The next section simulates the process outlined above. A total of 64
% sweeps are simulated and a range Doppler response is generated at the
% end.
%
% During the simulation, a spectrum analyzer is used to show the spectrum
% of each received sweep as well as its dechirped counterpart.

hspec = dsp.SpectrumAnalyzer('SampleRate',fs,...
    'PlotAsTwoSidedSpectrum',true,...
    'Title','Spectrum for received and dechirped signal',...
    'ShowLegend',true);
    

%%
% Next, run the simulation loop. 

%chirp个数
Nsweep = 64;
xr = complex(zeros(hwav.SampleRate*hwav.SweepTime,Nsweep));

for m = 1:Nsweep
    [radar_pos,radar_vel] = step(...
        hradarplatform,hwav.SweepTime);       % Radar moves during sweep
    [tgt_pos,tgt_vel] = step(hcarplatform,... 
        hwav.SweepTime);                      % Car moves during sweep
    x = step(hwav);                           % Generate the FMCW signal
    xt = step(htx,x);                         % Transmit the signal
    xt = step(hchannel,xt,radar_pos,tgt_pos,...
        radar_vel,tgt_vel);                   % Propagate the signal
    xt = step(hcar,xt);                       % Reflect the signal
    xt = step(hrx,xt);                        % Receive the signal
    xd = dechirp(xt,x);                       % Dechirp the signal
    
%    step(hspec,[xt xd]);                      % Visualize the spectrum
    
    %一共1100*64点数据
    xr(:,m) = xd;                             % Buffer the dechirped signal
end


%%
% From the spectrum scope, one can see that although the received signal is
% wideband (channel 1), sweeping through the entire bandwidth, the
% dechirped signal becomes narrowband (channel 2). 

%% Range and Doppler Estimation
% Before estimating the value of the range and Doppler, it may be a good
% idea to take a look at the zoomed range Doppler response of all
% 64 sweeps.

% hrdresp = phased.RangeDopplerResponse('PropagationSpeed',c,...
%     'DopplerOutput','Speed','OperatingFrequency',fc,'SampleRate',fs,...
%     'RangeMethod','FFT','SweepSlope',sweep_slope,...
%     'RangeFFTLengthSource','Property','RangeFFTLength',2048,...
%     'DopplerFFTLengthSource','Property','DopplerFFTLength',256);
% 
% clf;
% plotResponse(hrdresp,xr)                 % Plot range Doppler map
% axis([-v_max v_max 0 range_max])
% clim = caxis;

%%
% From the range Doppler response, one can see that the car in front is a
% bit more than 40 m away and appears almost static. This is expected
% because the radial speed of the car relative to the radar is only 4 km/h,
% which translates to a mere 1.11 m/s.
%
% There are many ways to estimate the range and speed of the target car.
% For example, one can choose almost any spectral analysis method to
% extract both the beat frequency and the Doppler shift. This example uses
% the root MUSIC algorithm to extract both the beat frequency and the
% Doppler shift.
%
% As a side note, although the received signal is sampled at 150 MHz so the
% system can achieve the required range resolution, after the dechirp, one
% only needs to sample it at a rate that corresponds to the maximum beat
% frequency. Since the maximum beat frequency is in general less than the
% required sweeping bandwidth, the signal can be decimated to alleviate the
% hardware cost. The following code snippet shows the decimation process.

% Dn = fix(fs/(2*fb_max));
% for m = size(xr,2):-1:1
%     xr_d(:,m) = decimate(xr(:,m),Dn,'FIR');
% end
% fs_d = fs/Dn;

%%
% To estimate the range, firstly, the beat frequency is estimated using the
% coherently integrated sweeps and then converted to the range.

% fb_rng = rootmusic(pulsint(xr_d,'coherent'),1,fs_d);
% rng_est = beat2range(fb_rng,sweep_slope,c)


%%
% Second, the Doppler shift is estimated across the sweeps at the range
% where the target is present.

% peak_loc = val2ind(rng_est,c/(fs_d*2));
% fd = -rootmusic(xr_d(peak_loc,:),1,1/tm);
% v_est = dop2speed(fd,lambda)/2


%%
% Note that both range and Doppler estimation are quite accurate.

%% Range Doppler Coupling Effect
% One issue associated with linear FM signals, such as an FMCW signal, is
% the range Doppler coupling effect. As discussed earlier, the target range
% corresponds to the beat frequency. Hence, an accurate range estimation
% depends on an accurate estimate of beat frequency. However, the presence
% of Doppler shift changes the beat frequency, resulting in a biased range
% estimation. 
%
% For the situation outlined in this example, the range error caused by the
% relative speed between the target and the radar is

% deltaR = rdcoupling(fd,sweep_slope,c)

%%
% This error is so small that we can safely ignore it. 
% 
% Even though the current design is achieving the desired performance, one
% parameter warrants further attention. In the current configuration, the
% sweep time is about 7 microseconds. Therefore, the system needs to sweep
% a 150 MHz band within a very short period. Such an automotive radar may
% not be able to meet the cost requirement. Besides, given the velocity of
% a car, there is no need to make measurements every 7 microseconds. Hence,
% automotive radars often use a longer sweep time. For example, the
% waveform used in [2] has the same parameters as the waveform designed in
% this example except a sweep time of 2 ms.
%
% A longer sweep time makes the range Doppler coupling more prominent. To
% see this effect, first reconfigure the waveform to use 2 ms as the sweep
% time.

% hwavtr = clone(hwav);
% release(hwavtr);
% tm = 2e-3;
% hwavtr.SweepTime = tm;
% sweep_slope = bw/tm;

%%
% Now calculate the range Doppler coupling. 

% deltaR = rdcoupling(fd,sweep_slope,c)

%%
% A range error of 1.14 m can no longer be ignored and needs to be
% compensated. Naturally, one may think of doing so following the same
% procedure outlined in earlier sections, estimating both range and
% Doppler, figuring out the range Doppler coupling from the Doppler shift,
% and then remove the error from the estimate.
%
% Unfortunately this process doesn't work very well with the long sweep
% time. The longer sweep time results in a lower sampling rate across the
% sweeps, thus reducing the radar's capability of unambiguously detecting
% high speed vehicles. For instance, using a sweep time of 2 ms, the
% maximum unambiguous speed the radar system can detect using the
% traditional Doppler processing is

% v_unambiguous = dop2speed(1/(2*tm),lambda)/2

%%
% The unambiguous speed is only 0.48 m/s, which means that the relative
% speed, 1.11 m/s, cannot be unambiguously detected. This means that not
% only the target car will appear slower in Doppler processing, the range
% Doppler coupling also cannot be correctly compensated.
%
% One way to resolve such ambiguity without Doppler processing is to adopt
% a triangle sweep pattern. Next section shows how the triangle sweep
% addresses the issue.




%速度距离解算
close all;
%中频信号
figure;
subplot(2,2,1);
plot(real(xt));
title('Receive the signal');

%dechirp之后的信号
subplot(2,2,2);
plot(real(xd));
title('Dechirp the signal');

sweep_number = 64;
samlpe_points = 2048;

fft1 = (zeros(samlpe_points,sweep_number));
fft2 = (zeros(samlpe_points,sweep_number));
%对实部做fft
for n=1:sweep_number
    fft1(:,n) = fft(real(xr(:,n)),samlpe_points);
    fft1_f = ((1:samlpe_points)-1)*150e6/samlpe_points;
    hold on;
end

for n=1:samlpe_points
    fft2(n,:) = abs(fft(fft1(n,:),sweep_number));
end

%距离维频谱图
freq = ((1:samlpe_points)-1)*fs/samlpe_points;
range = freq*c*tm/bw/2;
subplot(2,2,3)
plot(range(1:samlpe_points/2),abs(fft1(1:samlpe_points/2)));


%速度维频谱图
dopple_f = ((1:sweep_number)-1)*1/tm/sweep_number*lambda/2*3.6;
for n = 1:10   
    subplot(2,2,4)
    plot(dopple_f(1:sweep_number/2),abs(fft2(n,1:sweep_number/2)));
    hold on;
end


figure;
x = linspace(0,1/tm*lambda/2*3.6,sweep_number);
y = linspace(0,fs/sweep_slope*c/2,samlpe_points);


mesh(x,y,fft2);
ylim([0 range_max]);
xlim([0 330]);

 ylabel('distance:m');
 xlabel('velocity:Km/h');
 title('FMCW target detection')
