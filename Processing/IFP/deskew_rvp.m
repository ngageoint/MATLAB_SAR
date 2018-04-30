function [pulses_deskew, t_start] = deskew_rvp(pulses, sampling_rate, chirp_rate, pad, time_shift)
%DESKEW_RVP Perform range deskew on a set of pulses
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Range deskew is also known as residual video phase (RVP) compensation.
% It is essentionally just a linearly frequency-dependent time shift.  See
% appendix C of Carrara, Goodman, and Majewski for a more in-depth
% analysis.
%
% USAGE:
%   [pulses_deskew, t_start] = deskew_rvp(pulses, sampling_rate, chirp_rate, pad)
%
% INPUTS:
%   pulses        - required : pulses to be deskewed
%   sampling_rate - required : sampling rate for the system
%   chirp_rate    - required : chirp rate for the system
%   pad           - optional : pad data so that no data is shifted off the
%                   ends.  Otherwise output samples have exactly the same
%                   time indices as input samples.  (Default false)
%   time_shift    - optional : vector of per pulse values that represent a
%                   constant time shift for all frequencies within a pulse
%                   to accomodate "dejitter" processing.  (Default is zero
%                   for all pulses.)
%
% OUTPUTS:
%   pulses_deskew - required : deskewed pulses
%   t_start       - optional : starting time of first output sample,
%                   relative to the first input sample.  If PAD input
%                   parameter is true, then this will be non-zero.
%
% VERSION:
%   1.0
%     - Sean Hatch 20080222
%     - initial version
%   1.1
%     - Wade Schwartzkopf 20120329
%     - Added pad option
%   1.2
%     - Wade Schwartzkopf 20130805
%     - Added dejitter option
%
% TODO:

% Calculate the size of FFT and output data
[num_samples, num_pulses] = size(pulses);
total_time_shift = abs((sampling_rate^2)/chirp_rate); % Range of time shifts (in samples) from min to max frequency
num_fft = num_samples + total_time_shift; % Pad enough to assure that deskew does not wrap around
num_fft = 2^nextpow2(num_fft); % Bump up to fast FFT time
if exist('pad','var')&&pad
    pre_pad = ceil(total_time_shift/2); % Don't allow shift to wrap around before first sample
    output_samples = num_samples + ceil(total_time_shift);
else
    pre_pad = 0;
    output_samples = num_samples;
end

% Starting time of output data, with respect to input data.  We have likely
% shifted some frequencies before the start of the original data samples.
% If PAD is true, then we are returning the data shifted there, as well as
% samples at the original time indices.
t_start = - pre_pad/sampling_rate;

% Compute deskew phase
% Start by defining frequency axis.  num_fft steps of
% sampling_rate/num_fft.  Must be zero (no deskew shift) at DC.
N = num_fft; % Shorten notation in next line
f = ifftshift((sampling_rate/N)*(-ceil((N-1)/2):(((N-1)-ceil((N-1)/2))))).';
phase = -(f.^2)/chirp_rate; % Quadratic (deskew) term
phase = repmat(phase, 1, num_pulses);
if exist('time_shift','var')
    phase = phase + 2*f*time_shift(:).'; % Linear (dejitter) term
end
phase = exp(pi*1i*phase);

% Pad the pulse data appropriately
pulses = [repmat(complex(0), pre_pad, num_pulses); ...
          pulses; ...
          repmat(complex(0), (num_fft - num_samples) - pre_pad, num_pulses)];

% Apply the phase correction for deskew in FFT domain
pulses_deskew = ifft(fft(pulses) .* phase);

% Remove the zero padding
pulses_deskew = pulses_deskew(1:output_samples, :);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////