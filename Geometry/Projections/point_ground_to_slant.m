function outpoints = point_ground_to_slant(points, metadata, varargin)
% POINT_GROUND_TO_SLANT Projects focus (output) plane points to the
% range-Doppler slant plane points.  This is the inverse of the
% point_slant_to_ground function.
%
% This function is mostly just a wrapper to support legacy code that
% depended on an older version of complex image-to-ground projection.  It
% now just calls the SICD sensor model code.  This function varies from the
% main sensor model code in the toolbox (point_image_to_ground) in that it
% takes lat/long/hae coordinates as opposed to ECF, and it uses 1-based
% indexing into pixels (the MATLAB convention) as opposed to 0-based
% indexing (the SICD convention.)
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

outpoints = round(point_scene_to_image(geodetic_to_ecf(points), metadata) + 1);

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
