clc;clear all;close all;

fc = 77e9;%center frequency
c = 3e8;%light velocity in free space
lambda = c/fc;% wave length

target_dist = 50; 
target_speed = 110;
target_speed = target_speed*1000/3600;

target_dist1 = 125; 
target_speed1 = 60;
target_speed1 = target_speed1*1000/3600;
range_max = 200;%maximum of detection distance

tm = 10*2*range_max/c;%sweep time
range_rs = 1;% range resolution
bw = c/2/range_rs;%sweep bandwidth
sweep_slope = bw/tm;% slope

fr_max = sweep_slope*2*range_max/c;% beat frequency corresponding to the maximum range
v_max = 230*1000/3600; %maximum speed is about 230 km/h
fd_max = 2*v_max/lambda; % maximum doppler frequency
fb_max = fr_max+fd_max ; %maximum beat frequency
fs = 10*fb_max; %sampling frequency
t = 0:1/fs:tm; % sampling time 

Nsweep = 64; % pulses number
Range_FFT_N = 2048;
xr =  (zeros(length(t),Nsweep));
xr0 =  (zeros(length(t),Nsweep));
xr1 =  (zeros(length(t),Nsweep));
xr2 =  (zeros(length(t),Nsweep));
range_fft =  (zeros(Range_FFT_N,Nsweep));
speed_fft =  (zeros(Range_FFT_N,Nsweep));

for m = 1:1:Nsweep
    % ԭʼ�ź�
    xr0(:,m) = exp(1i*2*pi*fc*t+1i*pi*sweep_slope*t.^2);
    
%     spectrogram(xr0(:,m),128,120,128,1e6,'yaxis')
    % echo signal_1
    xr1(:,m) = exp(1i*2*pi*(    (2*target_dist/lambda+2*target_speed*(m-1)*tm/lambda)+...
                                      (    (2*target_dist/c*sweep_slope+2*target_speed/lambda+2*target_speed/c*sweep_slope*(m-1)*tm)*t)+...
                                           2*target_speed/c*sweep_slope*t.^2)...
                                      ); 
    % echo signal_2
    xr2(:,m) = exp(1i*2*pi*(  (  2*target_dist1/lambda+2*target_speed1*(m-1)*tm/lambda)+...
                                        (( 2*target_dist1/c*sweep_slope+2*target_speed1/lambda+2*target_speed1/c*sweep_slope*(m-1)*tm)*t)+...
                                           2*target_speed1/c*sweep_slope*t.^2)...
                                      ); 

    xr = xr1(:,m)+xr2(:,m);
%     figure(1)
%     subplot(211);plot(t,real(xr1(:,m)),t,real(xr2(:,m)))
%     subplot(212);plot(t,real(xr))
    range_fft(:,m) =fft(real(xr),2048);  
end

for n = 1:Range_FFT_N/2
    speed_fft(n,:) = abs(fft(((range_fft(n,:))),64));% velocity
end

%����άƵ��ͼ
figure;
range = ((1:2048)-1)*fs/2048*c*tm/bw/2;
for i=1:64
    plot(range,abs(range_fft(1:2048,i)));hold on;
end

%�ٶ�άƵ��ͼ
figure;
speed = ((1:64))*1/tm/64*lambda/2*3.6;
for i=1:Range_FFT_N/2
    plot(speed,abs(speed_fft(i,:)),'-b*');hold on;
end

figure;
mesh(speed,range,speed_fft);
ylim([0 250]);
xlim([0 330]);
view(3);
ylabel('distance:m');
xlabel('velocity:Km/h');
title('FMCW target detection');

