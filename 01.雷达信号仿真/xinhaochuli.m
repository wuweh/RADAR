%%��ʾ������Դ�ԡ��״��źŴ��������̡�

%%%%%%%%%%%%%%%%%%%%%%%% ȥб���������� %%%%%%%%%%%%%%%%%%%%%%%%% 
clc;clear all;close all;

B=10e6;%���� 10MHz
tp=10e-6;%���� 10us
k=B/tp;%LFM ϵ��
fs=50e6; R0=3e3;R1=2000;R2=3500;R=5000; c=3e8;
f0=60e6;
N=round(2*R/c*fs); fft_N=2^nextpow2(N); t=linspace(0,2*R/c,N);


%%%%%%%%%%%%%%%%%%%%%%% �ο��ź� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sref=exp(2i*pi*f0*t).*exp(1i*pi*k*t.^2);
Sref=exp(2i*pi*f0*t+1i*pi*k*t.^2);
figure
plot(real(Sref));
%%%%%%%%%%%%%%%%%%%%%%%%%%% �ز��ź� %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Sb0=exp(1j*pi*k*(t-2*R0/c).^2).*exp(2j*pi*f0*(t-2*R0/c)); 
Sb1=exp(1j*pi*k*(t-2*R1/c).^2).*exp(2j*pi*f0*(t-2*R1/c)); 
Sb2=exp(1j*pi*k*(t-2*R2/c).^2).*exp(2j*pi*f0*(t-2*R2/c)); 
Sb=Sb0+Sb1+Sb2;
figure
plot(real(Sb));   %������Ŀ����ɵĻز��ź�ʵ��

%%%%%%%%%%%%%%%%%%%%%%%%%%% ��Ƶ�ź� %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
SSb=Sref.*conj(Sb);         %ȥб��ʱ���źţ�����Ƶ����Ƶ
figure
plot(abs(SSb));         %abs(SSb)����adc����ֵ

spectrum=fft(SSb,fft_N);    %ȥб��Ƶ���ź� 
figure
plot(abs(spectrum));

f=fs*(0:fft_N-1)/fft_N-fs/2;    %��-fs/2 �� fs/2 
f=f*c*tp/2/B;               %˲ʱƵ�ʶ�Ӧ�ľ���
sf=exp(-j*pi/k*f.^2);       %�˲������亯�� 
SSb=spectrum.*sf;       %��Ƶ��ȥ����Ť����ʵ����ѹ����ȥ RVP figure;
SSb=fftshift(SSb); 
SSb1=ifft(SSb);         %�����˾���Ť���� RVP ��ʱ���ź� 
figure
subplot(211);
plot(f,db(abs(SSb)/max(SSb)))
xlabel('����/m');
grid on
subplot(212);
plot(f,abs(SSb))
xlabel('����/m');
grid on