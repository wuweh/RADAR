clc;
clear all;
close all;

fc = 77e9;  %77GHz
c = 3e8;      %����
lambda = c/fc;  %����

%����Զ������ȷ��ɨƵʱ��
range_max = 200;                                        %��Զ������
tm = 5.5*range2time(range_max,c);           %ɨƵʱ��

% ����ֱ��ʣ�range resolution) = c * Tp / 2 
% range resolution = C/(2*B);
% ���� B ���Ƶ��� 1/Tp, ���Դ���Խ�󣬷ֱ���Խ�ߡ�
range_res = 1;
bw = range2bw(range_res,c);    %�ɾ���ֱ���ȷ��ɨƵ����
sweep_slope = bw/tm;                %ɨƵб�� = ɨƵ����/ɨƵʱ��

% The function computes 2*range_max*sweep_slope/c 
%����Զ���������������ƵƵ��
fr_max = range2beat(range_max,sweep_slope,c);

%������ٶ�����Ķ�����Ƶ�ʣ�fd_max = 2*v_max/lambda;
%2*v_max����Ϊ�����������౳��ʻ
v_max = 70;
fd_max = speed2dop(2*v_max,lambda);

%�ɴ˿��Եó����յ�����ƵƵ��
fb_max = fr_max+fd_max;

%ȷ������Ƶ��
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


%����FMCW���β���
hwav = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw, 'SampleRate',fs);

% �γɲ���
s = step(hwav);
%ʵ�����źŵķ��ȣ�������
% subplot(211); plot(0:1/fs:tm-1/fs,real(s));
% xlabel('Time (s)'); ylabel('Amplitude (v)');
% title('FMCW signal'); axis tight;

%spectrogram: �ö�ʱ����Ҷ�任
%16��hanmming windows��
%0 sample overlap
% 4096��FFT
%����Ƶ��Ϊfs
%���䲨��Ƶ��ͼ
% subplot(212); spectrogram(s,32,16,4096,fs,'yaxis');
% title('FMCW signal spectrogram');

car_dist = 30;
car_speed = 10;

%Ŀ�공����RCSֵ
car_rcs = db2pow(min(10*log10(car_dist)+5,20));

%cָ���ǹ��٣�fcָ�����ز�Ƶ��
hcar = phased.RadarTarget('MeanRCS',car_rcs,'PropagationSpeed',c,...
    'OperatingFrequency',fc);

%����Ŀ�공���ĳ�ʼλ��[x,y,z]�ͳ�ʼ�ٶ�[vx,vy,vz]
hcarplatform = phased.Platform('InitialPosition',[car_dist;5;0.5],...
    'Velocity',[car_speed;0;0]);

% The propagation model is assumed to be free space.

hchannel = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

%t���߿׾���С m^2
ant_aperture = 6.06e-4;                         % in square meter
%ת��Ϊ��������
ant_gain = aperture2gain(ant_aperture,lambda);  % in dB

tx_ppower = db2pow(5)*1e-3;                     % in watts
tx_gain = 9+ant_gain;                           % in dB

rx_gain = 15+ant_gain;                          % in dB
rx_nf = 4.5;                                    % in dB

htx = phased.Transmitter('PeakPower',tx_ppower,'Gain',tx_gain);
hrx = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,...
    'SampleRate',fs);

%�����ٶ� 100km/h�ͳ�ʼλ����Ϣ
radar_speed = 20;
hradarplatform = phased.Platform('InitialPosition',[0;0;0.5],...
    'Velocity',[radar_speed;0;0]);

% hspec = dsp.SpectrumAnalyzer('SampleRate',fs,...
%     'PlotAsTwoSidedSpectrum',true,...
%     'Title','Spectrum for received and dechirped signal',...
%     'ShowLegend',true);

%chirp����
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

    %һ��1100*64������
    xr(:,m) = xd;                             % Buffer the dechirped signal
end


%�ٶȾ������
%��Ƶ�ź�
figure;
subplot(2,2,1);
plot(real(xt));
title('Receive the signal');

%dechirp֮����ź�
subplot(2,2,2);
plot(real(xd));
title('Dechirp the signal');

sweep_number = 64;
Range_FFT_P = 2048;

fft_R = (zeros(Range_FFT_P,sweep_number));
fft_V = (zeros(Range_FFT_P/2,sweep_number));

%��ʵ����fft
for n=1:sweep_number
    fft_R(:,n) = fft(real(xr(:,n)),Range_FFT_P);
    hold on;
end
%����άƵ��ͼ
freq = ((1:Range_FFT_P)-1)*fs/Range_FFT_P;
range = freq*c*tm/bw/2;
subplot(2,2,3)
plot(range(1:Range_FFT_P/2),abs(fft_R(1:Range_FFT_P/2)));

%�ٶ�άFFT
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
 
 
 
 %����Ŀ��
