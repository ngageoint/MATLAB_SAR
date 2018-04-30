function [ x, y, z ] = wgs_84_norm( x, y, z )
%WGS_84_NORM Computes the normal vector to the WGS_84 ellipsoid at a given
%point in ECEF space
%
% [ norm_xyz ] = wgs_84_norm( ecef_xyz )
%
% where ecef_xyz is a 3-by-n array of n ECEF vectors (XYZ, meters), and the
% resulting norm_xyz is of the same size.
%
% [ norm_x, norm_y, norm_z ] = wgs_84_norm( ecef_x, ecef_y, ecef_z )
%
% where the X, Y, and Z components are each passed separately.  ecef_x,
% ecef_y, and ecef_z must all be of the same size.  The resulting norm_x,
% norm_y, and norm_z will all be of the same size as well.
%
% Written by: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Define constants
a=6378137; % Semi-major (equatorial) axis of WGS_84 model
b=6356752.314245179; % Semi-minor (polar) axis of WGS_84 model

% Handle different forms in input arguments
if nargin==3 % Componentwise inputs, separate arguments for X,Y,Z
    % Nothing to do.  Processing uses this form.
elseif size(x,1)==3 % Array of 3-element vectors
    z = x(3,:);
    y = x(2,:);
    x = x(1,:);
elseif numel(x)==3 % Horizontal 3-element vector
    z = x(3);
    y = x(2);
    x = x(1);
else
    error('WGS_84_NORM:INVALID_INPUTS', 'Invalid inputs.');
end

% Calculate normal vector
x = x/(a^2);
y = y/(a^2);
z = z/(b^2);
% Make into unit vector
mag = sqrt(x.^2 + y.^2 + z.^2);
x = x./mag;
y = y./mag;
z = z./mag;

if nargout < 2 % Matrix/vector form, rather than componentwise form, was requested
    x = [x; y; z];
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////