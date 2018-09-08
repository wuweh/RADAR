%%%%%%%%%%%%%%%%%%%%%%%% 窄带加窗处理 %%%%%%%%%%%%%%%%%%%%%%%% 
clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%% 去斜加窗处理 %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
B=10e6;%带宽 10MHz 
tp=10e-6;%脉宽 10us 
u=B/tp;%LFM 系数 
fs=50e6;%fs>=2*B/tp*tau 
R0=3000;%初始距离 
R=4500;%距离波门
c=3e8; 
f0=60e6;%载频
N=round(2*R/c*fs); 
fft_N=2^nextpow2(N); 
t=linspace(0,2*R/c,N); 
f=fs*(0:fft_N-1)/fft_N-fs/2;%从-fs/2 到 fs/2

%%%%%%%%%%%%%%%%%%%%%%%%%% 参考信号 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Sref=exp(1i*pi*u*t.^2);
% figure
% plot(real(Sref));
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 回波信号 %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Sb=rectpuls(t-2*R0/c,tp).*exp(1j*pi*u*(t-2*R0/c).^2); 
% figure
% plot(real(Sb));
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 混频信号 %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ssb=Sref.*conj(Sb);
% figure
% plot(abs(ssb),'b');
% hold on

%% 时域加窗 
hamming=[zeros(749,1)',hamming(502)',zeros(249,1)']; 
ssb0=hamming.*ssb;
% figure
% plot(abs(ssb0),'r')

spectrum_ssb0=fft(ssb0,fft_N); %加窗后fft结果
spectrum_ssb=fft(ssb,fft_N);    %未加窗fft结果
f=f*c*tp/2/B;%瞬时频率对应的距离
figure; %%图 6 
subplot(2,1,1)
hold on
plot(f,db(abs(fftshift(spectrum_ssb))/max(fftshift(spectrum_ssb))),'b') 
plot(f,db(abs(fftshift(spectrum_ssb0))/max(fftshift(spectrum_ssb0))), 'r')

subplot(2,1,2)
hold on
plot(abs(fftshift(spectrum_ssb)),'b');
plot(abs(fftshift(spectrum_ssb0)),'r');
hold off