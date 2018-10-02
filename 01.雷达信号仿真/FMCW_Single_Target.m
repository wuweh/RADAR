clc;
clear all;
close all;

fc = 77e9;  %77GHz
c = 3e8;      %光速
lambda = c/fc;  %波长

%由最远检测距离确定扫频时间
range_max = 200;                                        %最远检测距离
tm = 5.5*range2time(range_max,c);           %扫频时间

% 距离分辨率（range resolution) = c * Tp / 2 
% range resolution = C/(2*B);
% 带宽 B 近似等于 1/Tp, 所以带宽越大，分辨率越高。
range_res = 1;
bw = range2bw(range_res,c);    %由距离分辨率确定扫频带宽
sweep_slope = bw/tm;                %扫频斜率 = 扫频带宽/扫频时间

% The function computes 2*range_max*sweep_slope/c 
%由最远检测距离带来的最大差频频率
fr_max = range2beat(range_max,sweep_slope,c);

%由最大速度引起的多普勒频率：fd_max = 2*v_max/lambda;
%2*v_max是因为两车可能是相背行驶
v_max = 70;
fd_max = speed2dop(2*v_max,lambda);

%由此可以得出最终的最大差频频率
fb_max = fr_max+fd_max;

%确定采样频率
fs = max(2*fb_max,bw);


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


%生成FMCW波形参数
hwav = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw, 'SampleRate',fs);

% 形成波形
s = step(hwav);
%实部即信号的幅度！！！！
% subplot(211); plot(0:1/fs:tm-1/fs,real(s));
% xlabel('Time (s)'); ylabel('Amplitude (v)');
% title('FMCW signal'); axis tight;

%spectrogram: 用短时傅里叶变换
%16阶hanmming windows；
%0 sample overlap
% 4096点FFT
%采样频率为fs
%发射波形频谱图
% subplot(212); spectrogram(s,32,16,4096,fs,'yaxis');
% title('FMCW signal spectrogram');

car_dist = 30;
car_speed = 10;

%目标车辆的RCS值
car_rcs = db2pow(min(10*log10(car_dist)+5,20));

%c指的是光速，fc指的是载波频率
hcar = phased.RadarTarget('MeanRCS',car_rcs,'PropagationSpeed',c,...
    'OperatingFrequency',fc);

%定义目标车辆的初始位置[x,y,z]和初始速度[vx,vy,vz]
hcarplatform = phased.Platform('InitialPosition',[car_dist;5;0.5],...
    'Velocity',[car_speed;0;0]);

% The propagation model is assumed to be free space.

hchannel = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

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

%本车速度 100km/h和初始位置信息
radar_speed = 20;
hradarplatform = phased.Platform('InitialPosition',[0;0;0.5],...
    'Velocity',[radar_speed;0;0]);

% hspec = dsp.SpectrumAnalyzer('SampleRate',fs,...
%     'PlotAsTwoSidedSpectrum',true,...
%     'Title','Spectrum for received and dechirped signal',...
%     'ShowLegend',true);

%chirp个数
Nsweep = 64;
xr = complex(zeros(hwav.SampleRate*hwav.SweepTime,Nsweep));

for m = 1:Nsweep
    [radar_pos,radar_vel] = step(hradarplatform,hwav.SweepTime);       % Radar moves during sweep
    [tgt_pos,tgt_vel] = step(hcarplatform,hwav.SweepTime);                      % Car moves during sweep
    x = step(hwav);                           % Generate the FMCW signal
    xt = step(htx,x);                         % Transmit the signal
    xt = step(hchannel,xt,radar_pos,tgt_pos, radar_vel,tgt_vel);                   % Propagate the signal
    xt = step(hcar,xt);                       % Reflect the signal
    xt = step(hrx,xt);                        % Receive the signal
    xd = dechirp(xt,x);                       % Dechirp the signal

    %    step(hspec,[xt xd]);                      % Visualize the spectrum

    %一共1100*64点数据
    xr(:,m) = xd;                             % Buffer the dechirped signal
end


%速度距离解算
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
Range_FFT_P = 2048;

fft_R = (zeros(Range_FFT_P,sweep_number));
fft_V = (zeros(Range_FFT_P/2,sweep_number));

%对实部做fft
for n=1:sweep_number
    fft_R(:,n) = fft(real(xr(:,n)),Range_FFT_P);
    hold on;
end
%距离维频谱图
freq = ((1:Range_FFT_P)-1)*fs/Range_FFT_P;
range = freq*c*tm/bw/2;
subplot(2,2,3)
plot(range(1:Range_FFT_P/2),abs(fft_R(1:Range_FFT_P/2)));

%速度维FFT
for n=1:1:Range_FFT_P/2
    fft_V(n,:) = fft(real(fft_R(n,:)),sweep_number);
    fft_V(n,:) = fftshift(fft_V(n,:));
end

dopple_f = ((1:sweep_number)-33)*1/tm/sweep_number*lambda/2;
subplot(2,2,4)
for n=1:Range_FFT_P/4
    plot(dopple_f(1:64),abs(fft_V(n,1:64)));
    hold on;
end

% figure;
% mesh(dopple_f(1:sweep_number/2),range(1:Range_FFT_P/2),abs(fft_V(1:Range_FFT_P/2,1:sweep_number/2)));
% view(3);
% ylabel('distance:m');
% xlabel('velocity:Km/h');
% title('FMCW target detection')
 
 
 
 %解算目标
