function [ boolean_out ] = is_normalized_sicd( sicd_meta_struct, dim)
%IS_NORMALIZED_SICD Test whether SAR data is already normalized in a given dimension
%
% Normalized complex data will have
%    1) Centered and deskewed frequency support
%    2) Uniform weighting
%    3) FFT sign of -1
%
% Author: Wade Schwartzkopf, NGA/R
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if ~exist('dim','var')
    dim = 1;
end

if dim==1 % slow-time
    if isfield(sicd_meta_struct,'Grid') && isfield(sicd_meta_struct.Grid,'Col')
        sicd_grid_struct = sicd_meta_struct.Grid.Col;
    else
        boolean_out = false; return;
    end
else % fast-time
    if isfield(sicd_meta_struct,'Grid') && isfield(sicd_meta_struct.Grid,'Row')
        sicd_grid_struct = sicd_meta_struct.Grid.Row;
    else
        boolean_out = false; return;
    end
end

% Several reasons that we might need to applied normalization
needs_deskew = isfield(sicd_grid_struct,'DeltaKCOAPoly') && ...
    any(sicd_grid_struct.DeltaKCOAPoly(:)~=0);
not_uniform_weight = (isfield(sicd_grid_struct,'WgtType') && ...
    isfield(sicd_grid_struct.WgtType,'WindowName') && ...
    ~strcmp(sicd_grid_struct.WgtType.WindowName,'UNIFORM')) || ...
    (isfield(sicd_grid_struct,'WgtFunct') && ...
    any(diff(sicd_grid_struct.WgtFunct)));
needs_fft_sgn_flip = isfield(sicd_grid_struct,'Sgn') && ...
    sicd_grid_struct.Sgn(1) == 1;

boolean_out = ~needs_deskew && ~not_uniform_weight && ~needs_fft_sgn_flip;

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////