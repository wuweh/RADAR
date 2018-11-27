clear all;
clc;
close all;

fc = 77e9;                      %77GHz
c = 3e8;                         %����
lambda = c/fc;              %����

% range2time = 2*r/c.
%����Զ������ȷ��ɨƵʱ��
range_max = 200;                                        %��Զ������
tm = 5.5*range2time(range_max,c);           %ɨƵʱ����������Ϊ5.5����Զ���봫��ʱ��

%ɨƵ����ȡ���ھ���ֱ���
% ����ֱ��ʣ�range resolution) = c * Tp / 2 
% ���� B ���Ƶ��� 1/Tp, ���Դ���Խ�󣬷ֱ���Խ�ߡ�
% range resolution = C/(2*B);
range_res = 1;
bw = range2bw(range_res,c);    %�ɾ���ֱ���ȷ��ɨƵ����
sweep_slope = bw/tm;           %ɨƵб�� = ɨƵ����/ɨƵʱ��


% fr_max =  2*range_max*sweep_slope/c 
%����Զ���������������ƵƵ��
fr_max = range2beat(range_max,sweep_slope,c);


%������ٶ�����Ķ�����Ƶ�ʣ�fd_max = 2*v_max/lambda;
%2*v_max����Ϊ�����������౳��ʻ
v_max = 230*1000/3600;
fd_max = speed2dop(2*v_max,lambda);

%�ɴ˿��Եó����յ�����ƵƵ��
fb_max = fr_max+fd_max;

%ȷ������Ƶ��
fs = max(2*fb_max,bw);

%77G�״�ϵͳ����
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

%д�벨�β�����ɨƵʱ�䣻ɨƵ����������
hwav = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw, 'SampleRate',fs);
plot(hwav);
% �γɲ���
s = step(hwav);
%ʵ�����źŵķ��ȣ�������
subplot(211); 
plot(0:1/fs:tm-1/fs,real(s));
xlabel('Time (s)'); ylabel('Amplitude (v)');
title('FMCW signal'); axis tight;

%spectrogram: �ö�ʱ����Ҷ�任
%32��hanmming windows��
%16 sample overlap
% 32��FFT
%����Ƶ��Ϊfs
%���䲨��Ƶ��ͼ
subplot(212); spectrogram(s,32,16,256,fs,'yaxis');
title('FMCW signal spectrogram');

%�״�ɢ�����(RCS)������Ծ�����ƣ�
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

%t���߿׾���С m^2
ant_aperture = 6.06e-4;                         % in square meter
%ת��Ϊ��������
ant_gain = aperture2gain(ant_aperture,lambda);  % in dB

%���书��
tx_ppower = db2pow(5)*1e-3;                     % in watts
%��������
tx_gain = 9+ant_gain;                           % in dB

%��������
rx_gain = 15+ant_gain;                          % in dB
rx_nf = 4.5;                                    % in dB

htx = phased.Transmitter('PeakPower',tx_ppower,'Gain',tx_gain);
hrx = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,'SampleRate',fs);

%�����ٶ� 100km/h�ͳ�ʼλ����Ϣ
radar_speed = 0*1000/3600;
hradarplatform = phased.Platform('InitialPosition',[0;0;0.5],'Velocity',[radar_speed;0;0]);

hspec = dsp.SpectrumAnalyzer('SampleRate',fs,'PlotAsTwoSidedSpectrum',true,'Title',...
'Spectrum for received and dechirped signal','ShowLegend',true);

% Next, run the simulation loop. 

%chirp����
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
    xt = xt_1 +xt_2;   %�����Ļز��ź�
    %���ݽ������ߵ����棬�õ������ź�
    xt = step(hrx,xt);                        % Receive the signal
    %�����Ե�Ƶ
    xd = dechirp(xt,x);                       % Dechirp the signal
    %     step(hspec,[xt xd]);                      % Visualize the spectrum
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
Speed_FFT_P = 128;

fft_R = (zeros(Range_FFT_P,sweep_number));
fft_V = (zeros(Range_FFT_P/2,sweep_number));

%��ʵ����fft
for n=1:sweep_number
    fft_R(:,n) = fft(real(xr(:,n)),Range_FFT_P);
    hold on;
end
%����άƵ��ͼ
freq = ((1:Range_FFT_P)-1)*fs/Range_FFT_P;
range = freq(1:Range_FFT_P/2)*c*tm/bw/2;
subplot(2,2,3)
plot(range(1:Range_FFT_P/2),abs(fft_R(1:Range_FFT_P/2)));

%�ٶ�άFFT
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
