clc;clear all;close all
%% parameters setting
 
B = 135e6;                       % Sweep bandwidth
T = 36.5e-6;                     % Sweep time
N = 512;                         % Sample length
L = 128;                         % Chirp total
c = 3e8;                         % Speed of light
f0 = 76.5e9;                     % Start frequency
NumRFFT = 512;                   % Range fft length
NumDFFT = 128;                   % Dopler fft length
rangeRes = c/2/B;                % Range resolution
velRes = c/2/f0/T/NumDFFT;       % Velocity resolution
maxRange = rangeRes*NumRFFT/2;   % Max range
maxVel = velRes*NumDFFT/2;       % Max velocity
tarR = [50; 90];                 % Target range
tarV = [3; 20];                  % Target velocity

fs = N/T;
t = 0:1/fs:T;
k = B/T;

figure;
xr0(:,1) = exp(1i*pi*k*(t.^2));

subplot(211); plot(real(xr0(:,1)));
xlabel('Time'); ylabel('Amplitude (v)');
title('FMCW signal'); axis tight;

% spectrogram: 用短时傅里叶变换
% 16阶hanmming windows；
% 0 sample overlap
% 32点FFT
% 采样频率为fs
% 发射波形频谱图
subplot(212); spectrogram(xr0(:,1),16,8,512,N/T,'yaxis');
title('FMCW signal spectrogram');



%% generate receive signal
S1 = zeros(L,N);
for l = 1:L
    for n = 1:N
        S1(l,n) = exp(1i*2*pi*...
                                        (...
                                           ( 2*B*(tarR(1)+tarV(1)*T*l)/(c*T) + (2*f0*tarV(1))/c )*T/N*n...
                                       +    (2*f0*(tarR(1)+tarV(1)*T*l)/c)...
                                   )...
                               );
    end
end
 
S2 = zeros(L,N);
for l = 1:L
    for n = 1:N
        S2(l,n) = exp(1i*2*pi*...
                                        (...
                                           ( 2*B*(tarR(2)+tarV(2)*T*l)/(c*T)+(2*f0*tarV(2))/c )*T/N*n...
                                       +    (2*f0*(tarR(2)+tarV(2)*T*l)/c)...
                                       )...
                             );
    end
end
 
sigReceive = S1+S2;
 
%% range win processing
sigRWin = zeros(L,N);
for ii = 1:L
    sigRWin(ii,:) = sigReceive(ii,:).*hamming(N)';
end
 
figure;
%% range fft processing
sigRfft = zeros(L,NumRFFT);
range = ((1:NumRFFT)-1)*(N/T)/NumRFFT*c*T/B/2;
for ii = 1:L
    sigRfft(ii,:) = fft(sigRWin(ii,:),NumRFFT);
    subplot(121)
    plot(range,abs(sigRfft(ii,:)));hold on
end
 
%% doppler win processing
sigDWin = zeros(L,NumRFFT);
for ii = 1:NumRFFT
    sigDWin(:,ii) = sigRfft(:,ii).*hamming(L);
end
 
%% doppler fft processing
sigDfft = zeros(NumDFFT,NumRFFT);
speed = ((1:NumDFFT))*1/T/NumDFFT*(c/f0)/2;
for ii = 1:NumRFFT
    sigDfft(:,ii) = fft(sigDWin(:,ii),NumDFFT);
    subplot(122)
    plot(speed,abs(sigDfft(:,ii)));hold on
end
 
figure
sigDfftShift = fftshift(sigDfft(:,1:NumRFFT/2),1);
sigDfftShift_mag = abs(sigDfftShift);
mesh(sigDfftShift_mag);
