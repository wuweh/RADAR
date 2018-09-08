function [xr,xr_unmixed] = helperFMCWSimulate(Nsweep,waveform,...
    radarmotion,carmotion,transmitter,channel,cartarget,receiver)
% This function helperFMCWSimulate is only in support of FMCWExample. It
% may be removed in a future release.

%   RSWEEP =
%   helperFMCWSimulate(NSWEEP,WAVEFORM,RADARMOTION,CARMOTION,TRANSMITTER,
%   CHANNEL,CARTARGET,RECEIVER) returns the simulated sweep train RSWEEP. 
%
%   The input parameters are:
%       NSWEEP:             number of sweeps
%       WAVEFORM:           waveform object
%       RADARMOTION:        platform object for the radar
%       CARMOTION:          platform object for target car
%       TRANSMITTER:        transmitter object
%       CHANNEL:            propagation channel object
%       CARTARGET:          target car object
%       RECEIVER:           receiver object
%
%   The rows of RSWEEP represent fast time and its columns represent slow
%   time (pulses). When the pulse transmitter uses staggered PRFs, the
%   length of the fast time sequences is determined by the highest PRF.

%   Copyright 2010-2016 The MathWorks, Inc.

release(waveform);
release(radarmotion);
release(transmitter);
release(receiver);
release(carmotion);
release(channel);
release(cartarget);

if isa(waveform,'phased.MFSKWaveform')
    sweeptime = waveform.StepTime*waveform.StepsPerSweep;
else
    sweeptime = waveform.SweepTime;
end
Nsamp = round(waveform.SampleRate*sweeptime);

xr = complex(zeros(Nsamp,Nsweep));
xr_unmixed = xr;

Ntgt = numel(cartarget.MeanRCS);
for m = 1:Nsweep
    % Update radar and target positions
    [radar_pos,radar_vel] = radarmotion(sweeptime);
    [tgt_pos,tgt_vel] = carmotion(sweeptime);

    % Transmit FMCW waveform
    sig = waveform();
    txsig = transmitter(sig);

    % Propagate the signal and reflect off the target
    rxsig = complex(zeros(Nsamp,Ntgt));
    for n = 1:Ntgt
        rxsig(:,n) = channel(txsig,radar_pos,tgt_pos(:,n),radar_vel,tgt_vel(:,n));
    end
    rxsig = cartarget(rxsig);
    
    % Dechirp the received radar return
    rxsig = receiver(sum(rxsig,2));
    xd = dechirp(rxsig,sig);
    xr_unmixed(:,m) = rxsig;
    xr(:,m) = xd;
end
