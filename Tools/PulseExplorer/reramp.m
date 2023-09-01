function [ rfPulse dt ] = reramp( pulses, sampleRate, derampRate)
%RERAMP Upsample and reramp the input complex pulse signal
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% USAGE:
%   [ rfPulse t ] = reramp( pulses, sampleRate, derampRate );
%
% INPUTS:
%   pulses     - matrix of pulse data ordered column-wise
%   sampleRate - sample rate (Hz)
%   derampRate - deramp rate (Hz/s)
%
% OUTPUTS:
%   rfPulse - optional: upsampled and reramped pulse data
%   dt      - optional: corresponding sample timing for the 'rfPulse' data
%
% VERSION:
%   1.0 20071210 Ralph Fiedler
%     - initial version
%   1.1 20111211 Wade Schwartzkopf
%     - streamlined input parameters


% Compute the required upsample factor to hold reramped data without
% aliasing
% - The sample timing is defined as the inverse of the sample rate, e.g.,
%     actualDeltaT = 1 / sampleRate
% - The desired sample timing is defined by the inverse of the chirp
%   bandwidth plus the instanteous bandwidth from the sampling rate, e.g.,
%     desiredDeltaT = 1 / BW = 1 / (( derampRate * rcvWindowLength ) +
%     sampleRate)
% - The upsample factor is just the ratio of these two delta times, e.g.,
%     upsampleFactor = actualDeltaT / desiredDeltaT = BW / sampleRate
nsamp = size(pulses,1);
rcv_window_length = nsamp / sampleRate;
deramp_bw = abs(rcv_window_length * derampRate);
upsampleFactor = ( deramp_bw + sampleRate ) / sampleRate;

% "Oversample" factor brings you beyond what is strictly required for this
% data.  Here we only bring to an integer sample size, but we could also
% include a factor to bring up to a fast FFT size.
n = round(upsampleFactor * nsamp);
overSampleFactor = n / upsampleFactor / nsamp;

% Zero-mean the pulse data
pulses = bsxfun(@minus, pulses, mean(pulses));

% Use a FFT-based upsample method
rfPulse = interpft( pulses, n );

% Compute the sample time interval
dt = 1 / sampleRate / upsampleFactor / overSampleFactor;

% Compute the sample times for the reramp
% Set the origin, t=0, to the middle of the receive window
t = (( 0 : n-1 )' * dt) - (rcv_window_length / 2);

% Compute the reramp function
% This will cause the center of receive window (in time) to be centered (in
% frequency) at DC:
ramp = exp( pi * 1i * derampRate * t .* t );
% This will shift in frequency so that at the center of the receive window
% (in time), it will be centered (in frequency) at the center of the new
% upsampled bandwidth (effectively performing an FFTSHIFT):
% ramp = exp( pi * 1i * ( (sampleRate * upsampleFactor * overSampleFactor * t) + (derampRate * t .* t) ) );

% Apply the reramp
rfPulse = bsxfun(@times, rfPulse, ramp);

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////