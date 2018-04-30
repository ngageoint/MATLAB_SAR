function [ DeltaKCOAPoly, az_coords_m, rg_coords_m, fft_sgn ] = deskewparams( sicd_meta, dim )
%DESKEWPARAMS Extract from SICD structure the parameters for deskewmem.m
%    [ DeltaKCOAPoly, az_coords_m, rg_coords_m ] = deskewparams( sicd_meta, dim )
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% DeltaKCOA polynomial
if dim==1 && isfield(sicd_meta,'Grid') && isfield(sicd_meta.Grid,'Col') &&...
        isfield(sicd_meta.Grid.Col,'DeltaKCOAPoly') % deskew in azimuth
    DeltaKCOAPoly = sicd_meta.Grid.Col.DeltaKCOAPoly;
elseif dim==2 && isfield(sicd_meta,'Grid') && isfield(sicd_meta.Grid,'Row') &&...
        isfield(sicd_meta.Grid.Row,'DeltaKCOAPoly') % deskew in range
    DeltaKCOAPoly = sicd_meta.Grid.Row.DeltaKCOAPoly;
else % DeltaKCOA is optional.  If not filled out, assume to be zero.
    DeltaKCOAPoly = 0;
end
% Vectors describing range and azimuth distances from SCP (in meters) for rows and columns
az_coords_m = (double(0:(sicd_meta.ImageData.NumCols-1)) + ...
    double(sicd_meta.ImageData.FirstCol) - double(sicd_meta.ImageData.SCPPixel.Col)) * ...
    sicd_meta.Grid.Col.SS;
rg_coords_m = (double(0:(sicd_meta.ImageData.NumRows-1)) + ...
    double(sicd_meta.ImageData.FirstRow) - double(sicd_meta.ImageData.SCPPixel.Row)) * ...
    sicd_meta.Grid.Row.SS;
% FFT sign required to transform data to spatial frequency domain
if dim==1 && isfield(sicd_meta,'Grid') && isfield(sicd_meta.Grid,'Col') &&...
        isfield(sicd_meta.Grid.Col,'Sgn') % azimuth
    fft_sgn = sicd_meta.Grid.Col.Sgn;
elseif dim==2 && isfield(sicd_meta,'Grid') && isfield(sicd_meta.Grid,'Row') &&...
        isfield(sicd_meta.Grid.Row,'Sgn') % range
    fft_sgn = sicd_meta.Grid.Row.Sgn;
else
    fft_sgn = -1; % Default, since this is more common
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////