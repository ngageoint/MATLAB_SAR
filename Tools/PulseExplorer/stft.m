function [ yy, t, f ] = stft( data, sample_rate, center )
%STFT Short Time Fourier Transform
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% This routine generates and (optionally) plots the stft of a data series.
% This is done by doing the fft on a sliding window through the data.
%
% Requires Signal Processing Toolbox.
%
% INPUTS:
%   data        - required: complex time domain data
%   sample_rate - required: sample rate of data (Hz)
%   center      - optional: specify for FFT shift about the center
% 
% OUTPUTS:
%   See spectrogram documentation
 
% define stft window parameters
nfft = 4 * ( 2^nextpow2( sqrt(length(data)) ) );
window = kaiser( floor(0.97*nfft), 5 );
noverlap = floor( 0.9 * nfft );

% compute stft
[yy,f,t]=spectrogram(data, window, noverlap, nfft, sample_rate);
% Older syntax:
% [output_data,f,t]=specgram(data, nfft, sample_rate, window, noverlap);
yy = abs(yy);
if center&&~isreal(data)
    yy = fftshift(yy,1);
end
% Another way to do this, that looks better in some cases
% [yy2,f2,t2]=fsst(data,sample_rate,'yaxis');

if nargout<1
    figure;
    imagesc(t,f,yy);
    axis xy;
    title( 'Spectrogram' );
    xlabel( 'Fast-Time (s)' );
    ylabel( 'Frequency' );
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////