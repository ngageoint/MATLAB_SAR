function [x, y, z] = geodetic_to_ecf(lat, lon, alt)
%GEODETIC_TO_ECF Convert geodetic coordinates to ECF
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
%
% Convert geodetic latitude, longitude, and altitude to ECF (Earth 
% Centered Fixed) coordinates.
%
% USAGE:
%   pos_ecf = geodetic_to_ecf(pos_lla)
%   [pos_ecf_x, pos_ecf_y, pos_ecf_z] = geodetic_to_ecf(lat, lon, alt)
%
% INPUTS:
%   pos_lla - required : geodetic latitude, longitude, and altitude   [deg, deg, m]
%
% OUTPUTS:
%   pos_ecf - required : ecf x, y, z coordinates                      [m, m, m]
%
% NOTES:
%   Zhu, J. Conversion of Earth-centered, Earth-fixed coordinates to 
%   geodetic coordinates. IEEE Transactions on Aerospace and Electronic
%   Systems, 30, 3 (July 1994), 957-962.
%
% VERSION:
%   1.0
%     - Sean Hatch 20070911
%     - initial version
%   1.1
%     - Wade Schwartzkopf 20130708
%     - vectorized and componentwise data handling
%
% TODO:

% define constants
e2 = 6.6943799901377997e-3; % eccentricity squared of Earth (WGS 84 value)
a = 6378137.0;              % semimajor radius of the Earth (WGS 84 value)

% Handle different forms in input arguments
if nargin==3 % Componentwise inputs, separate arguments for lat,lon,alt
    % Nothing to do.  Processing uses this form.
elseif size(lat,1)==3 % Array of 3-element vectors
    alt = lat(3,:);
    lon = lat(2,:);
    lat = lat(1,:);
elseif numel(lat)==3 % Horizontal 3-element vector
    alt = lat(3);
    lon = lat(2);
    lat = lat(1);
else
    error('WGS_84_NORM:INVALID_INPUTS', 'Invalid inputs.');
end

% calculate distance to surface of ellipsoid
R = a ./ sqrt(1.0 - e2 .* sind(lat) .* sind(lat));

% calculate coordinates
x = (R + alt) .* cosd(lat) .* cosd(lon);
y = (R + alt) .* cosd(lat) .* sind(lon);
z = (R + alt - e2 .* R) .* sind(lat);

if nargout < 2 % Matrix/vector form, rather than componentwise form, was requested
    x = [x; y; z];
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////