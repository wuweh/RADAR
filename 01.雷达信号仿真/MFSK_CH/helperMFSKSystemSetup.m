function [fmcwwaveform,target,tgtmotion,channel,transmitter,receiver,...
    sensormotion,c,fc,lambda,fs,fr_max] = helperMFSKSystemSetup
% This function helperMFSKSystemSetup is only in support of MFSKExample. It
% may be removed in a future release.

%   Copyright 2014-2016 The MathWorks, Inc.

% System parameter
fc = 77e9;      % operating frequency
c = 3e8;        % propagation speed
lambda = c/fc;  % wavelength

tm = 0.001;     % sweep time

range_res = 1;              % range resolution
bw = range2bw(range_res,c); % bandwidth
sweep_slope = bw/tm;        % sweep slope

range_max = 200;
fr_max = range2beat(range_max,sweep_slope,c);

v_max = 230*1000/3600;
fd_max = speed2dop(2*v_max,lambda);

fb_max = fr_max+fd_max;


fs = max(2*fb_max,bw);


fmcwwaveform = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw,...
    'SampleRate',fs,'SweepDirection','Triangle');

car_dist = 50;
car_speed = 20*1000/3600;
car_rcs = db2pow(min(10*log10(car_dist)+5,20));

truck_dist = 55;
truck_speed = 40*1000/3600;
truck_rcs = db2pow(min(10*log10(truck_dist)+5,20));

bike_dist = 65;
bike_speed = -60*1000/3600;
bike_rcs = db2pow(min(10*log10(bike_dist)+5,20));

tgtpos = [[car_dist;0;0],[truck_dist;0;0],[bike_dist;0;0]];
tgtvel = [[car_speed;0;0],[truck_speed;0;0],[bike_speed;0;0]];
tgtmotion = phased.Platform('InitialPosition',tgtpos,'Velocity',tgtvel);

tgtrcs = [car_rcs,truck_rcs,bike_rcs];
target = phased.RadarTarget('MeanRCS',tgtrcs,'PropagationSpeed',c,...
    'OperatingFrequency',fc);

channel = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fc,'SampleRate',fs,'TwoWayPropagation',true);

ant_aperture = 6.06e-4;                         % in square meter
ant_gain = aperture2gain(ant_aperture,lambda);  % in dB

tx_ppower = db2pow(5)*1e-3;                     % in watts
tx_gain = 9+ant_gain;                           % in dB

rx_gain = 15+ant_gain;                          % in dB
rx_nf = 4.5;                                    % in dB

transmitter = phased.Transmitter('PeakPower',tx_ppower,'Gain',tx_gain);
receiver = phased.ReceiverPreamp('Gain',rx_gain,'NoiseFigure',rx_nf,...
    'SampleRate',fs);

radar_speed = 60*1000/3600;
sensormotion = phased.Platform('Velocity',[radar_speed;0;0]);


