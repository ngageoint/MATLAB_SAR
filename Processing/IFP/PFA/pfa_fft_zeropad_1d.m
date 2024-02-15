function [data] = pfa_fft_zeropad_1d(data, sample_rate)
%PFA_FFT_ZEROPAD_1D Performs the zeropadding and FFT required by pfa_mem
%and pfa_file.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
    
% Setup zeropad
zeropad = zeros(floor(size(data,1)*sample_rate), size(data,2), class(data));
start_index = floor(size(data,1)*(sample_rate-1)/2);

zeropad(start_index + (1:size(data,1)),:) = data; % Insert data into zeropad
zeropad = ifftshift(zeropad, 1); % Move DC (center point) to index 1
data = fftshift(fft(zeropad,[],1),1); % Actual FFT

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
