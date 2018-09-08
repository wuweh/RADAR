%基于数学原理的雷达信号仿真

clc;
clear all;
close all;

%target_1
target_dist = 50; %distance  should be less than 200
target_speed = 80;
target_speed = target_speed*1000/3600;% speed should be less than 230

%target_2
target_dist1 = 125; %distance  should be less than 200
target_speed1 = 100;
target_speed1 = target_speed1*1000/3600;% speed should be less than 230

% radar parameter
fc = 77e9;%center frequency
c = 3e8;%light velocity in free space
lambda = c/fc;% wave length

range_max = 200;%maximum of detection distance
tm = 10*2*range_max/c%sweep time

range_rs = 0.5;% range resolution
bw = c/2/range_rs%sweep bandwidth
sweep_slope = bw/tm;% slope

fr_max = sweep_slope*2*range_max/c;% beat frequency corresponding to the maximum range

v_max = 230*1000/3600; %maximum speed is about 230 km/h
fd_max = 2*v_max/lambda; % maximum doppler frequency

fb_max = fr_max+fd_max  %maximum beat frequency

% fs = 10*fb_max%sampling frequency
fs = 20e6;
t = 0:1/fs:tm; % sampling points
FFT_R_N = 512;
FFT_V_N = 64;
f0 = 77e9-0.5*bw;

Nsweep = 64; % pulses number
for m = 1:1:Nsweep
    xr1(:,m) = exp(1i*2*pi*((2*target_dist/lambda+2*target_speed*(m-1)*tm/lambda)+...
    ((2*target_dist/c*sweep_slope+2*target_speed/lambda+2*target_speed/c*sweep_slope*(m-1)*tm)*t)+...
    2*target_speed/c*sweep_slope*t.^2));% echo signal_1

    xr2(:,m) = exp(1i*2*pi*((2*target_dist1/lambda+2*target_speed1*(m-1)*tm/lambda)+...
    ((2*target_dist1/c*sweep_slope+2*target_speed1/lambda+2*target_speed1/c*sweep_slope*(m-1)*tm)*t)+...
    2*target_speed1/c*sweep_slope*t.^2));% echo signal_2

%     xr3(:,m) = exp(1i*(2*pi*f0*t+pi*sweep_slope*(t-(m-1)*tm).^2));
xr3(:,m) = exp(1i*(2*pi*f0*t));

    xr = xr1(:,m)+xr2(:,m);
    fft1(:,m) =fft(xr,FFT_R_N); % distance fft
end

figure;
subplot(121)
plot(real(xr3(:,m)))
subplot(122)
plot(abs(fft(xr3(:,m),512)))


figure;
subplot(221)
plot(real(xr1(1:256)))
title('中频时域信号')

for n = 1:1:FFT_R_N
    fft2(n,:) = abs(fft((fft1(n,:)),FFT_V_N));% velocity fft
end 

%距离维频谱图
subplot(222)
freq = ((1:FFT_R_N)-1)*fs/FFT_R_N;
range = freq*c*tm/bw/2;
plot(range,abs(fft1));
title('距离维FFT')

%速度维频谱图
subplot(223)
d_f = ((1:64)-1)*1/tm/64*lambda/2*3.6;
for n = 1:FFT_R_N
    plot(abs(fft2(n,:)));
    hold on;
end
title('速度维FFT')

subplot(224)
x = linspace(0,1/tm*lambda/2*3.6,Nsweep);
y = linspace(0,fs/sweep_slope*c/2,FFT_R_N);

mesh(x,y,fft2);
ylim([0 range_max]);
xlim([0 330]);
view(3);
ylabel('distance:m');
xlabel('velocity:Km/h');
title('FMCW target detection')