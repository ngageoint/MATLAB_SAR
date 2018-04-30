function totalRCS = RCS_Compute(complex_data, oversample_ratio, cal_sf, mask)
%RCS_COMPUTE Demonstrates area-based computation of calibrated RCS for an ROI
%
% USAGE:
%   totalRCS = RCS_Compute(complex_data, oversample_ratio, cal_sf, mask)
%
% INPUTS:
%   complex_data     - required : Complex valued SAR dataset in the image
%                                 domain.  First dimension is azimuth,
%                                 second range.  Third dimension could be
%                                 for multi-channel (i.e. polarimetric)
%                                 data.  This code assumes the frequency
%                                 support is centered and constant across
%                                 the image.
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
%
% Author: Tim Cox, NRL; Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Default parameter values
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

% Image domain computation of total RCS (Reference Adam Bryant, NGA/PL)
totalRCS = (1/prod(oversample_ratio)) * sum(sum(cal_sf.*abs(filtimg).^2));
totalRCS = squeeze(totalRCS); % Multi-channel (polarimetric) data
% The oversample ratio (or zeropad) factor is the ratio between the sum of
% squared (power detected) samples of an ideal sinc function (which is what
% we are measuring, at least for an ROI which we assume contains nearly all
% the energy of that sinc) and the peak of that ideal sinc^2 function
% (which is how the calibration constant is defined).
% We put the scale factor inside the sum so that it can be applied per
% pixel, rather than only as a constant.

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////