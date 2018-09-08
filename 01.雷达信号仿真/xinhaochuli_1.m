%%%%%%%%%%%%%%%%%%%%%%%% խ���Ӵ����� %%%%%%%%%%%%%%%%%%%%%%%% 
clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%% ȥб�Ӵ����� %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
B=10e6;%���� 10MHz 
tp=10e-6;%���� 10us 
u=B/tp;%LFM ϵ�� 
fs=50e6;%fs>=2*B/tp*tau 
R0=3000;%��ʼ���� 
R=4500;%���벨��
c=3e8; 
f0=60e6;%��Ƶ
N=round(2*R/c*fs); 
fft_N=2^nextpow2(N); 
t=linspace(0,2*R/c,N); 
f=fs*(0:fft_N-1)/fft_N-fs/2;%��-fs/2 �� fs/2

%%%%%%%%%%%%%%%%%%%%%%%%%% �ο��ź� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Sref=exp(1i*pi*u*t.^2);
% figure
% plot(real(Sref));
%%%%%%%%%%%%%%%%%%%%%%%%%%%% �ز��ź� %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Sb=rectpuls(t-2*R0/c,tp).*exp(1j*pi*u*(t-2*R0/c).^2); 
% figure
% plot(real(Sb));
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ��Ƶ�ź� %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ssb=Sref.*conj(Sb);
% figure
% plot(abs(ssb),'b');
% hold on

%% ʱ��Ӵ� 
hamming=[zeros(749,1)',hamming(502)',zeros(249,1)']; 
ssb0=hamming.*ssb;
% figure
% plot(abs(ssb0),'r')

spectrum_ssb0=fft(ssb0,fft_N); %�Ӵ���fft���
spectrum_ssb=fft(ssb,fft_N);    %δ�Ӵ�fft���
f=f*c*tp/2/B;%˲ʱƵ�ʶ�Ӧ�ľ���
figure; %%ͼ 6 
subplot(2,1,1)
hold on
plot(f,db(abs(fftshift(spectrum_ssb))/max(fftshift(spectrum_ssb))),'b') 
plot(f,db(abs(fftshift(spectrum_ssb0))/max(fftshift(spectrum_ssb0))), 'r')

subplot(2,1,2)
hold on
plot(abs(fftshift(spectrum_ssb)),'b');
plot(abs(fftshift(spectrum_ssb0)),'r');
hold off