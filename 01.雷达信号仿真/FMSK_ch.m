clc;clear all;close all;

A = 5;
f0 = 24e9;
fstep = 150e3;
tstep = 10e-6;
BW = 150e6;
N = 128;
T = N*tstep;
fs=100e6;
t = 1/fs:1/fs:tstep;

for n=1:N
    x1(:,n)= A*exp(-1i*2*pi*(f0+(n-1)*fstep)*t);
end

figure;
plot(real(x1(:,5)))
figure
a = fft(real(x1(:,47)),1024);
plot(abs(a(1:512)))