function range_profile = RCS_Range_Profile( complex_data, ...
    oversample_ratio, cal_sf, mask)
%RCS_RANGE_PROFILE Computes calibrated range profile
%
% USAGE:
%   [range_profile] = RCS_Compute(complex_data, oversample_ratio, cal_sf, mask)
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
%   range_profile    - optional : Range profile.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Apply shape mask to rectangular data
filtimg = complex_data.*repmat(mask,[1 1 size(complex_data,3)]);
if ~isscalar(cal_sf)
    % Multi-channel (polarimetric)
    cal_sf = repmat(cal_sf,[1 1 size(complex_data,3)]);
end

range_profile = (1/prod(oversample_ratio)) * sum(cal_sf .* abs(filtimg).^2,1);
if isvector(range_profile)
    range_profile = range_profile(:); % Assure in first dimension
else
    range_profile = squeeze(range_profile); % For multi-channel (polarimetric) data
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////