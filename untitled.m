%% FIR Highpass Filter
% Load |chirp.mat|. The file contains a signal, |y|, that has most of its
% power above |Fs/4|, or half the Nyquist frequency. The sample rate is
% 8192 Hz.

% Copyright 2015 The MathWorks, Inc.


%%
% Design a 34th-order FIR highpass filter to attenuate the components of
% the signal below |Fs/4|. Use a cutoff frequency of 0.48 and a Chebyshev
% window with 30 dB of ripple.
close all
clc

load chirp

t = (0:length(y)-1)/Fs;

bhi = fir1(34,0.48,'high',chebwin(35,30));
freqz(bhi,1)

%%
% Filter the signal. Display the original and highpass-filtered signals.
% Use the same _y_-axis scale for both plots.

outhi = filter(bhi,1,y);

for i=1:length(y)
    if i<=35
        outhi2(i)  = sum(bhi(1:i).*y(1:i)');
    else
        outhi2(i)  = sum(bhi(1:35).*y(i-35+1:i)');
    end
end

figure;
subplot(3,1,1)
plot(t,y)
title('Original Signal')
ys = ylim;

subplot(3,1,2)
plot(t,outhi)
title('Highpass Filtered Signal')
xlabel('Time (s)')
ylim(ys)

subplot(3,1,3)
plot(t,outhi2)
title('Highpass Filtered Signal')
xlabel('Time (s)')
ylim(ys)

%%
% Design a lowpass filter with the same specifications. Filter the signal
% and compare the result to the original. Use the same _y_-axis scale for
% both plots.

blo = fir1(34,0.48,chebwin(35,30));
figure;freqz(blo,1)

outlo = filter(blo,1,y);
figure;
subplot(2,1,1)
plot(t,y)
title('Original Signal')
ys = ylim;

subplot(2,1,2)
plot(t,outlo)
title('Lowpass Filtered Signal')
xlabel('Time (s)')
ylim(ys)