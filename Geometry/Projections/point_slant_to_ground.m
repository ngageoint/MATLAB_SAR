function pos = point_slant_to_ground(points, metadata, varargin)
% POINT_SLANT_TO_GROUND Projects range-Doppler slant plane points to 
% focus (output) plane.
%
% This function is mostly just a wrapper to support legacy code that
% depended on an older version of complex image-to-ground projection.  It
% now just calls the SICD sensor model code.  This function varies from the
% main sensor model code in the toolbox (point_image_to_ground) in that it
% returns lat/long/hae coordinates as opposed to ECF, and it uses 1-based
% indexing into pixels (the MATLAB convention) as opposed to 0-based
% indexing (the SICD convention.)
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

try
    pos = ecf_to_geodetic(point_image_to_ground(points - 1, metadata, varargin{:}));
catch
    pos = [];
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////