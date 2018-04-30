function [ profile_slow, profile_fast ] = RCS_ST_FT( ...
    complex_data, lookdir, oversample_ratio, cal_sf, mask, meta)
%RCS_ST_FT Computes total calibrated RCS for an ROI and slow/fast time RCS data profiles
%
% USAGE:
%   [totalRCS, profile_slow, profile_fast] = RCS_ST_FT(...
%      complex_data, lookdir, oversample_ratio, cal_sf, mask)
%
% INPUTS:
%   complex_data     - required : Complex valued SAR dataset in the image
%                                 domain.  First dimension is azimuth,
%                                 second range.  Third dimension could be
%                                 for multi-channel (i.e. polarimetric)
%                                 data.  This code assumes the frequency
%                                 support is centered and constant across
%                                 the image.
%   lookdir          - required : "Left" or "Right".
%   oversample_ratio - optional : Overample or zeropad factor. (Required
%                                 for calibrated RCS.)
%   cal_sf           - optional : Calibration scale factor (linear).
%                                 Either a constant or an array the same
%                                 size as complex_data with per-pixel
%                                 values. (Required for calibrated RCS.)
%   mask             - optional : Binary image which is ones over the
%                                 region of interest. (Default is an image
%                                 of all ones.)
%
% OUTPUTS:
%   totalRCS         - required : Total calibrated RCS of ROI (calibrated
%                                    if oversample_ratio and cal_sf were
%                                    given.)
%   profile_slow     - optional : Slow-time profile.
%   profile_fast     - optional : Fast-time profile.
%
% Author: Tim Cox, NRL; Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Default parameter values
if ~exist('lookdir','var')
    lookdir = 'R';
end
if ~exist('oversample_ratio','var')
    oversample_ratio = [1 1]; % RCS values will be uncalibrated if this is unknown
end
if ~exist('cal_sf','var')
    cal_sf = 1; % RCS values will be uncalibrated if this is unknown
end
if ~exist('mask','var')
    mask = ones(size(complex_data));
end

% Apply shape mask to rectangular data
filtimg = complex_data.*repmat(mask,[1 1 size(complex_data,3)]);
if ~isscalar(cal_sf)
    % Only works if the radiometric scale factors are same for both channels
    cal_sf = repmat(cal_sf,[1 1 size(complex_data,3)]);
end

% Inverse polar formatting breaks on some data, so it is disabled until we
% can clean some things up.
inverse_polar = 0;

% Slow/Fast-time profiles
% Determine fft size to use based on size of imagery
data_size = size(complex_data); % Multi-channel (polarimetric) data will have more than 2 dimensions
fftsize = 2.^nextpow2(data_size([1 2]))*4.0;
if inverse_polar % Inverse polar formatting is the most precise way to do this
    totalRCS = RCS_Compute(complex_data, oversample_ratio, cal_sf, mask);
    % TODO: Handle multi-channel (polarimetric) datasets
    % TODO: Check automatically for whether data is suitable for inverse polar formatting
    % TODO: Pass k_a and k_r back out so that slow time axes can be computed with them
    % TODO: Show inverse polar format grid in figure if asked
    ss = [meta.Grid.Col.SS,meta.Grid.Row.SS];
    imp_resp_bw = [meta.Grid.Col.ImpRespBW,meta.Grid.Row.ImpRespBW];
    [angle_rf_domain, k_a, k_r] = pfa_inv_mem(filtimg, meta.Grid.Row.KCtr, ss, imp_resp_bw, fftsize, -1);
    profile_slow = sum(abs(angle_rf_domain).^2,2); % Sum over constant angle
    profile_fast = sum(abs(angle_rf_domain).^2,1); % Sum over constant RF frequency
    % Its not clear on how to compute calibrated RCS with inverse polar
    % formatting in the same way we do with the polar format approximation,
    % because of the polar inscription and changed sample rate, so we just
    % set the mean to be equal to the image domain computed RCS.
    profile_slow = profile_slow * totalRCS / mean(profile_slow);
    profile_fast = profile_fast * totalRCS / mean(profile_fast);
else % Otherwise use the polar format approximation, where columns approximate time/azimuth angle and rows approximate receive frequency
    % Direct translation of CASE Executive rcstrans.f code
    nz_data_points = fix(fftsize./oversample_ratio); % Determine non-zero data region
    % cal_sf must go in fft, rather than profile computation as a single scalar, since it may be a per-pixel array
    rgcomp = fftshift(fft(filtimg.*sqrt(cal_sf),fftsize(1),1),1)/fftsize(1);
    azcomp = fftshift(fft(filtimg.*sqrt(cal_sf),fftsize(2),2),2)/fftsize(2);
    profile_slow = sum(abs(rgcomp).^2,2) .* nz_data_points(1).^2 ./oversample_ratio(2);
    profile_fast = sum(abs(azcomp).^2,1) .* nz_data_points(2).^2 ./oversample_ratio(1);
end

% The following is another numerically identical way to compute RCS but
% might be easier to understand:
% Disperse in azimuth, so that columns approximate time/azimuth angle (ignoring polar formatting)
% rgcomp = fftshift(fft(filtimg,fftsize(1),1),1)/sqrt(fftsize(1));
% Disperse in range, so that rows approximate receive frequency (ignoring polar formatting)
% azcomp = fftshift(fft(filtimg,fftsize(2),2),2)/sqrt(fftsize(2));
% For MATLAB's FFT (and FFTW which MATLAB uses), normalizing the FFT by the
% square root of fftsize maintains constant total power, so that:
% norm(filtimg(:))==norm(rgcomp(:))==norm(azcomp(:))
% Note that this normalization factor is specific to this FFT
% implementation.  For instance, IDL's FFT is returned normalized and
% requires no scale factor.
% profile_slow = sum(abs(rgcomp).^2,2) ... % Total power (uncalibrated) per time/azimuth bin
%     * cal_sf ... % Scale by calibration constant
%     * nz_data_points(1) ... % Energy has now been split across nz_data_points(1) time/azimuth bins
%     / prod(oversample_ratio); % See comments in image domain RCS section on oversample factor
% profile_fast = sum(abs(azcomp).^2,1) ... % Total power (uncalibrated) per frequency bin
%     * cal_sf ... % Scale by calibration constant
%     * nz_data_points(2) ... % Energy has now been split across nz_data_points(2) frequency bins
%     / prod(oversample_ratio); % See comments in image domain RCS section on oversample factor
% We could check to make sure that our computation in the range- and/or
% azimuth-compressed domains match the computation from the image domain.
% Average RCS in range/az-compressed should equal total RCS in image
% domain.  There may be some small difference if energy leaks into the
% zeropad areas.
% nz_data_points = nz_data_points + mod(nz_data_points,2); % force to be even
% guard_band = (fftsize - nz_data_points) / 2;
% nz_data_start = guard_band + 1;
% nz_data_stop  = guard_band + nz_data_points;
% error = (mean(profile_slow(nz_data_start(1):nz_data_stop(1)))-totalRCS)/totalRCS;

if upper(lookdir(1))=='R'
    profile_slow = profile_slow(end:-1:1,:,:); % Flip slow time data
end
profile_slow = squeeze(profile_slow); % For multi-channel (polarimetric) data that comes in a 3rd dimension

if isvector(profile_fast)
    profile_fast = profile_fast(:); % Assure in first dimension
else
    profile_fast = squeeze(profile_fast); % For multi-channel (polarimetric) data
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////