clear all;
clc;
close all;

fc = 77e9;                      %77GHz
c = 3e8;                         %光速
lambda = c/fc;              %波长

% range2time = 2*r/c.
%由最远检测距离确定扫频时间
range_max = 200;                                        %最远检测距离
tm = 5.5*range2time(range_max,c);           %扫频时间这里设置为5.5倍最远距离传播时间

%扫频带宽取决于距离分辨率
% 距离分辨率（range resolution) = c * Tp / 2 
% 带宽 B 近似等于 1/Tp, 所以带宽越大，分辨率越高。
% range resolution = C/(2*B);
range_res = 1;
bw = range2bw(range_res,c);    %由距离分辨率确定扫频带宽
sweep_slope = bw/tm;           %扫频斜率 = 扫频带宽/扫频时间


% fr_max =  2*range_max*sweep_slope/c 
%由最远检测距离带来的最大差频频率
fr_max = range2beat(range_max,sweep_slope,c);


%由最大速度引起的多普勒频率：fd_max = 2*v_max/lambda;
%2*v_max是因为两车可能是相背行驶
v_max = 230*1000/3600;
fd_max = speed2dop(2*v_max,lambda);

%由此可以得出最终的最大差频频率
fb_max = fr_max+fd_max;

%确定采样频率
fs = max(2*fb_max,bw);

%77G雷达系统参数
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

%写入波形参数：扫频时间；扫频带宽；采样率
hwav = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw, 'SampleRate',fs);
plot(hwav);
% 形成波形
s = step(hwav);
%实部即信号的幅度！！！！
subplot(211); 
plot(0:1/fs:tm-1/fs,real(s));
xlabel('Time (s)'); ylabel('Amplitude (v)');
title('FMCW signal'); axis tight;

%spectrogram: 用短时傅里叶变换
%32阶hanmming windows；
%16 sample overlap
% 32点FFT
%采样频率为fs
%发射波形频谱图
subplot(212); spectrogram(s,32,16,256,fs,'yaxis');
title('FMCW signal spectrogram');

%雷达散射截面(RCS)可由相对距离估计；
car_dist_1 = 50;                                %50m
car_speed_1 = 60*1000/3600;         %120km/h
car_rcs_1 = db2pow(min(10*log10(car_dist_1)+5,20));
hcar_1 = phased.RadarTarget('MeanRCS',car_rcs_1,'PropagationSpeed',c,'OperatingFrequency',fc);
hcarplatform_1 = phased.Platform('InitialPosition',[car_dist_1;0;0.5],'Velocity',[car_speed_1;0;0]);

car_dist_2 = 80;                                    % 20m
car_speed_2 = 180*1000/3600;             %50km/h
car_rcs_2 = db2pow(min(10*log10(car_dist_2)+5,20));
hcar_2 = phased.RadarTarget('MeanRCS',car_rcs_2,'PropagationSpeed',c,'OperatingFrequency',fc);
hcarplatform_2 = phased.Platform('InitialPosition',[car_dist_2;0;0.5],'Velocity',[car_speed_2;0;0]);

hchannel = phased.FreeSpace('PropagationSpeed',c,'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

%t天线孔径大小 m^2
ant_aperture = 6.06e-4;                         % in square meter
%转换为天线增益
ant_gain = aperture2gain(ant_aperture,lambda);  % in dB

%发射功率
tx_ppower = db2pow(5)*1e-3;                     % in watts
%发射增益
tx_gain = 9+ant_gain;                           % in dB

%接收增益
rx_gain = 15+ant_gain;                          % in dB
rx_nf = 4.5;                                    % in dB

htx = phased.Transmitter('PeakPower',tx_ppower,'Gain',tx_gain);
hrx = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,'SampleRate',fs);

%本车速度 100km/h和初始位置信息
radar_speed = 0*1000/3600;
hradarplatform = phased.Platform('InitialPosition',[0;0;0.5],'Velocity',[radar_speed;0;0]);

hspec = dsp.SpectrumAnalyzer('SampleRate',fs,'PlotAsTwoSidedSpectrum',true,'Title',...
'Spectrum for received and dechirped signal','ShowLegend',true);

% Next, run the simulation loop. 

%chirp个数
Nsweep = 64;
xr = complex(zeros(hwav.SampleRate*hwav.SweepTime,Nsweep));

for m = 1:Nsweep
    [radar_pos,radar_vel] = step(hradarplatform,hwav.SweepTime);                       % Radar moves during sweep
    [tgt_pos_1,tgt_vel_1] = step(hcarplatform_1, hwav.SweepTime);                      % Car moves during sweep
    [tgt_pos_2,tgt_vel_2] = step(hcarplatform_2, hwav.SweepTime);                      % Car moves during sweep
    x = step(hwav);                           % Generate the FMCW signal
    xt = step(htx,x);                         % Transmit the signal
    xt_1 = step(hchannel,xt,radar_pos,tgt_pos_1,radar_vel,tgt_vel_1);                   % Propagate the signal_1
    xt_1 = step(hcar_1,xt_1);                       % Reflect the signal_1
    xt_2 = step(hchannel,xt,radar_pos,tgt_pos_2,radar_vel,tgt_vel_2);                   % Propagate the signal_2
    xt_2 = step(hcar_2,xt_2);                      % Reflect the signal_2
    xt = xt_1 +xt_2;   %完整的回波信号
    %根据接受天线的增益，得到接收信号
    xt = step(hrx,xt);                        % Receive the signal
    %解线性调频
    xd = dechirp(xt,x);                       % Dechirp the signal
    %     step(hspec,[xt xd]);                      % Visualize the spectrum
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
Speed_FFT_P = 128;

fft_R = (zeros(Range_FFT_P,sweep_number));
fft_V = (zeros(Range_FFT_P/2,sweep_number));

%对实部做fft
for n=1:sweep_number
    fft_R(:,n) = fft(real(xr(:,n)),Range_FFT_P);
    hold on;
end
%距离维频谱图
freq = ((1:Range_FFT_P)-1)*fs/Range_FFT_P;
range = freq(1:Range_FFT_P/2)*c*tm/bw/2;
subplot(2,2,3)
plot(range(1:Range_FFT_P/2),abs(fft_R(1:Range_FFT_P/2)));

%速度维FFT
for n=1:1:Range_FFT_P/2
    fft_V(n,1:Speed_FFT_P) = fft((fft_R(n,:)),Speed_FFT_P);
    fft_V(n,:) = fftshift(fft_V(n,:));
end

dopple_f = ((1:Speed_FFT_P)-Speed_FFT_P/2-1)*1/Speed_FFT_P/tm*lambda/2*3.6;
subplot(2,2,4)
for n=1:Range_FFT_P/4
    plot(dopple_f(1:Speed_FFT_P),abs(fft_V(n,1:Speed_FFT_P)));
    hold on;
end

figure;
mesh(dopple_f,range(1:Range_FFT_P/2),abs(fft_V(1:Range_FFT_P/2,1:Speed_FFT_P)));
view(3);
ylabel('distance:m');
xlabel('velocity:Km/h');
title('FMCW target detection')
