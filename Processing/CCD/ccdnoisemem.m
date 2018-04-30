function [ ccd_out, phase_out ] = ccdnoisemem( reference_image, match_image, ...
    corr_window_size, reference_noise_var, match_noise_var )
%CCDNOISEMEM Coherent change detection
%    [ccd_out, phase_out] = ccdnoisemem(reference_image, match_image, corr_window_size, reference_noise_var, match_noise_var)
%
% Compares two images by performing coherent change detection.  As opposed
% to the traditional CCD equation, this version is a maximum likelihood
% estimator that considers the noise level in the data.  (The traditional
% CCD equation assumed zero noise.) This code implements CCD equation as
% described in Wahl, Jakowatz, Yocky, "A New Optimal Change Estimator for
% SAR Coherent Change Detection".  Assumes images (and all processing) fit
% into memory and that input images are already registered.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if isscalar(corr_window_size)
    corr_window_size = corr_window_size * [1 1];
end
if ~exist('reference_noise_var','var')
    reference_noise_var = 0; % Resort to traditional CCD
end
if ~exist('match_noise_var','var')
    match_noise_var = 0; % Resort to traditional CCD
end

conjf_times_g = conv2(conj(reference_image).*match_image, ones(corr_window_size), 'same');
f_squared = conv2(abs(reference_image).^2, ones(corr_window_size), 'same');
g_squared = conv2(abs(match_image).^2, ones(corr_window_size), 'same');
ccd_out = 2 * abs(conjf_times_g)./ ... % New CCD formulation
    (f_squared + g_squared - prod(corr_window_size) * (reference_noise_var + match_noise_var));
ccd_out(~isfinite(ccd_out)) = 0; % Fix divide by zeros
ccd_out(ccd_out<0 | ccd_out>1) = 1; % Since we are using 1 to denote uncertainty
if nargout>1
    phase_out = angle(conjf_times_g);
    phase_out(~isfinite(phase_out))=0;
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////