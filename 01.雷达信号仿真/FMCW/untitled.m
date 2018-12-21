%--------------------------------------------------------------------------
% This Script is edited by RadarSay
% Created : 2018/08/22
% Description : generate fmcw receive signal && range doppler processing
%--------------------------------------------------------------------------
 
clc
clear
 
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
 
%% generate receive signal
 
S1 = zeros(L,N);
for l = 1:L
    for n = 1:N
        S1(l,n) = exp(1i*2*pi*(((2*B*(tarR(1)+tarV(1)*T*l)/(c*T)+(2*f0*tarV(1))/c)*T/N*n+((2*f0)*(tarR(1)+tarV(1)*T*l))/c)));
    end
end
 
S2 = zeros(L,N);
for l = 1:L
    for n = 1:N
        S2(l,n) = exp(1i*2*pi*(((2*B*(tarR(2)+tarV(2)*T*l)/(c*T)+(2*f0*tarV(2))/c)*T/N*n+((2*f0)*(tarR(2)+tarV(2)*T*l))/c)));
    end
end
 
sigReceive = S1+S2;
 
%% range win processing
 
sigRWin = zeros(L,N);
for ii = 1:L
    sigRWin(ii,:) = sigReceive(ii,:).*hamming(N)';
end
 
%% range fft processing
 
sigRfft = zeros(L,NumRFFT);
for ii = 1:L
    sigRfft(ii,:) = fft(sigRWin(ii,:),NumRFFT);
end
 
%% doppler win processing
sigDWin = zeros(L,NumRFFT);
for ii = 1:NumRFFT
    sigDWin(:,ii) = sigRfft(:,ii).*hamming(L);
end
 
%% doppler fft processing
 
sigDfft = zeros(NumDFFT,NumRFFT);
for ii = 1:NumRFFT
    sigDfft(:,ii) = fft(sigDWin(:,ii),NumDFFT);
end
 
sigDfftShift = fftshift(sigDfft(:,1:NumRFFT/2),1);
sigDfftShift_mag = abs(sigDfftShift);
 mesh(sigDfftShift_mag);
image(sigDfftShift_mag)
surf(sigDfftShift_mag)