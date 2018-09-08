function [xr,xr_unmixed] = helperFMCWTwoRaySimulate(Nsweep,waveform,...
    radarmotion,carmotion,transmitter,txchannel,rxchannel,cartarget,receiver)
% This function helperFMCWTwoRaySimulate is only in support of FMCWExample.
% It may be removed in a future release.

%   RSWEEP =
%   helperFMCWTwoRaySimulate(NSWEEP,WAVEFORM,RADARMOTION,CARMOTION,
%   TRANSMITTER, TXCHANNEL,RXCHANNEL,CARTARGET,RECEIVER) returns the
%   simulated sweep train RSWEEP.
%
%   The input parameters are:
%       NSWEEP:             number of sweeps
%       WAVEFORM:           waveform object
%       RADARMOTION:        platform object for the radar
%       CARMOTION:          platform object for target car
%       TRANSMITTER:        transmitter object
%       TXCHANNEL:          propagation channel object from radar to car
%       RXCHANNEL:          propagation channel object from car to radar
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
release(txchannel);
release(rxchannel);
release(cartarget);

if isa(waveform,'phased.MFSKWaveform')
    sweeptime = waveform.StepTime*waveform.StepsPerSweep;
else
    sweeptime = waveform.SweepTime;
end
Nsamp = round(waveform.SampleRate*sweeptime);

xr = complex(zeros(Nsamp,Nsweep));
xr_unmixed = xr;

for m = 1:Nsweep

    % Update radar and target positions
    [radar_pos,radar_vel] = radarmotion(sweeptime);
    [tgt_pos,tgt_vel] = carmotion(sweeptime);

    % Transmit FMCW waveform
    sig = waveform();
    txsig = transmitter(sig);

    % Propagate the signal and reflect off the target
    txsig = txchannel(txsig,radar_pos,tgt_pos,radar_vel,tgt_vel);               % propagate the signal
    rxsig = cartarget(txsig);
    rxsig = rxchannel(rxsig,radar_pos,tgt_pos,radar_vel,tgt_vel);               % propagate the signal
    
    % Dechirp the received radar return
    rxsig = receiver(rxsig);
    xd = dechirp(rxsig,sig);
    xr_unmixed(:,m) = rxsig;
    xr(:,m) = xd;    
end
