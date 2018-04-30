function [ ccd_out, phase_out ] = ccdmem( reference_image, match_image, corr_window_size )
%CCDMEM Coherent change detection
%    ccd_out = ccdmem(reference_image, match_image, corr_window_size)
%
% Compares two images by performing coherent change detection.  Implements
% CCD equation as described in Jakowatz, et al., "Spotlight-mode Synthetic
% Aperture Radar: A Signal Processing Approach".  Assumes images (and all
% processing) fit into memory and that input images are already registered.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

conjf_times_g = conv2(conj(reference_image).*match_image, ones(corr_window_size), 'same');
f_squared = conv2(abs(reference_image).^2, ones(corr_window_size), 'same');
g_squared = conv2(abs(match_image).^2, ones(corr_window_size), 'same');
ccd_out = abs(conjf_times_g)./sqrt(f_squared.*g_squared); % Jakowatz, eq. 5.102
ccd_out(~isfinite(ccd_out))=0; % Fix divide by zeros
if nargout>1
    phase_out = angle(conjf_times_g);
    phase_out(~isfinite(phase_out))=0; % Fix divide by zeros
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////