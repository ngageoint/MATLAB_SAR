function [ weight_fun ] = estimate_weighting_mem( data, dim, zeropad )
%ESTIMATE_WEIGHTING_MEM Estimate the weighting applied to complex SAR data
%
%    weight_fun = estimate_weighting_mem(data, dim, zeropad)
%
%       Parameter name    Description
% 
%       data              Complex SAR data in an array.  Assumes that the
%                            frequency support of data is centered.
%       dim               Dimension over which to estimate weighting.
%       zeropad           Zeropad factor in the dimension is which the 
%                            weighting is being estimated. (Default = 1)
%       weight_fun        Function handle that generates weighting.  Takes
%                            a single input parameter, which is the number
%                            of elements in the resulting weighting vector.
%
% Author: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if ~exist('dim','var')
    dim = 1;
end
if ~exist('zeropad','var')
    zeropad = 1;
end
other_dim = 3-dim;

%% Estimate weighting
weighting_estimate = sum(abs(fft(data,[],dim)),other_dim);
weighting_estimate = weighting_estimate(:); % Assure vertical vector
trim = floor((numel(weighting_estimate)/zeropad)/2);
zeropad_trim = weighting_estimate([1:trim (end-trim):end]);
% Smooth data, by fitting a raised cosine
raised_cos_fit = fft(zeropad_trim);
fft_elems = 3; % Enough for hamming, hanning, blackman (won't work for uniform though)
raised_cos_fit(fft_elems:(end-fft_elems+2)) = 0;
raised_cos_fit = ifft(abs(raised_cos_fit)); % The abs() centers the weighting
% Another way to smooth data:
% filter_size=ceil(numel(zeropad_trim)/16); % Wide averaging filter of arbitrary length
% filter_size=filter_size+mod(filter_size+1,2); % Make odd
% weighting_smoothed = conv(padarray(zeropad_trim,floor(filter_size/2),'circular'),...
%     ones(filter_size,1)/filter_size,'valid');
linear_fit = polyval(polyfit((1:numel(zeropad_trim))',zeropad_trim,1),...
    (1:numel(zeropad_trim))');
if norm(raised_cos_fit-zeropad_trim) < ... % Is raised cosine more likely
        norm(linear_fit-zeropad_trim)      % than linear fit?
    raised_cos_fit = raised_cos_fit/max(raised_cos_fit); % Make max 1
    weight_fun = @(x) fftshift(interpft(raised_cos_fit,x));
else
    weight_fun = []; % Data that matches linear fit is likely uniformly weighted
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////