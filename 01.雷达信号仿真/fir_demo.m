clc;clear all;close all;
load chirp
t = (0:length(y)-1)/Fs;
bhi = fir1(34,0.48,'high',chebwin(35,30));
freqz(bhi,1)

%%
% Filter the signal. Display the original and highpass-filtered signals.
% Use the same _y_-axis scale for both plots.

outhi = filter(bhi,1,y);

subplot(2,1,1)
plot(t,y)
title('Original Signal')
ys = ylim;

subplot(2,1,2)
plot(t,outhi)
title('Highpass Filtered Signal')
xlabel('Time (s)')
ylim(ys)

%%
% Design a lowpass filter with the same specifications. Filter the signal
% and compare the result to the original. Use the same _y_-axis scale for
% both plots.

blo = fir1(34,0.48,chebwin(35,30));

outlo = filter(blo,1,y);

subplot(2,1,1)
plot(t,y)
title('Original Signal')
ys = ylim;

subplot(2,1,2)
plot(t,outlo)
title('Lowpass Filtered Signal')
xlabel('Time (s)')
ylim(ys)