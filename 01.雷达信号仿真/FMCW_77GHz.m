%基于数学原理的雷达信号仿真(未完成)

clc;
clear all;
close all;

FFT_R_N = 512;
FFT_V_N = 64;
f0 = 77e9;

% radar parameter
fc = 77e9;%center frequency
c = 3e8;%light velocity in free space
lambda = c/fc;% wave length

range_max = 200;%maximum of detection distance
tm = 10*2*range_max/c;%sweep time

range_rs = 0.5;% range resolution
bw = c/2/range_rs;%sweep bandwidth
sweep_slope = bw/tm;% slope

fr_max = sweep_slope*2*range_max/c;% beat frequency corresponding to the maximum range

v_max = 230*1000/3600; %maximum speed is about 230 km/h
fd_max = 2*v_max/lambda; % maximum doppler frequency

fb_max = fr_max+fd_max;  %maximum beat frequency

% fs = 10*fb_max%sampling frequency
fs = 20e6;
t = 1:1:FFT_R_N; % sampling points

Nsweep = 64; % pulses number

fs = 500e6;
A = 10;
n = 0;

%target1
%target_1
target_dist = 50; %distance  should be less than 200
target_speed = 80;
target_speed = target_speed*1000/3600;% speed should be less than 230
delt_t1 = 2*target_dist/c;

%target_2
target_dist1 = 125; %distance  should be less than 200
target_speed1 = 100;
target_speed1 = target_speed1*1000/3600;% speed should be less than 230

bw1 = 0;
fm_1 = 2*bw*target_dist/c/tm;
fa_1 = 2*pi*(2*target_dist/c*(f0+bw/2)-2*bw*target_dist/tm/c/c);
for i=0:1/fs:tm
    n=n+1;
    s(n) = A*cos(2*pi*((f0+bw/2)*i + 1/2/tm*bw*i^2));
    s1(n) = A*cos(2*pi*((f0+bw/2)*(i-delt_t1) + 1/2/tm*bw*((i-delt_t1) )^2));
    delt_s(n) = A*cos(2*pi*((fm_1+bw1/2)*i + 1/2/tm*bw1*i^2)+fa_1);
%     delt_s(n) = s(n) - s1(n);
end
figure
plot(s);
figure;
plot(delt_s);

n = 0;
for i=1:100:6600
    n = n + 1;
    delt_s1(n) = delt_s(i);
end
figure;
plot(delt_s1);

range_fft = fft(delt_s1,128);
figure
plot(abs(range_fft))







