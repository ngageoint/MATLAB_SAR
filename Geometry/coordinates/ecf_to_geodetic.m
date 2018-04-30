function [lat, lon, alt] = ecf_to_geodetic(x, y, z)
%ECF_TO_GEODETIC Convert ECF coordinates to geodetic
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Convert ECF (Earth Centered Fixed) coordinates to geodetic latitude, 
% longitude, and altitude.
%
% USAGE:
%   pos_lla = ecf_to_geodetic(pos_ecf)
%   [lat, lon, alt] = geodetic_to_ecf(pos_ecf_x, pos_ecf_y, pos_ecf_z)
%
% INPUTS:
%   pos_ecf - required : ecf x, y, z coordinates                      [m, m, m]
%
% OUTPUTS:
%   pos_lla - required : geodetic latitude, longitude, and altitude   [deg, deg, m]
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
e2 = 6.6943799901377997e-3;  % eccentricity squared of Earth (WGS 84 value)
a = 6378137.0;               % semimajor radius of the Earth (WGS 84 value)

% calculate derived constants
e4 = e2 .* e2;
ome2 = 1.0 - e2;
a2 = a .* a;
b = a .* sqrt(ome2);         % semiminor radius of the Earth
b2 = b .* b;
e_b2 = (a2 - b2) ./ b2;

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

% calculate intermediates
z2 = z .* z;
r2 = (x .* x) + (y .* y);
r = sqrt(r2);

% Check for invalid solution
valid = ((a .* r) .* (a .* r) + (b .* z) .* (b .* z) > (a2 - b2) .* (a2 - b2));
lon = nan(size(x)); % Default values for invalid solutions
lat = nan(size(x));
alt = nan(size(x));

% calculate longitude
lon(valid) = atan2(y(valid), x(valid))*180/pi; % atan2d not available until MATLAB 2012b

% calculate intermediates
F = 54.0 .* b2 .* z2;
G = r2 + ome2 .* z2 - e2 .* (a2 - b2);
c = e4 .* F .* r2 ./ (G .* G .* G);
s = (1.0 + c + sqrt(c .* c + 2 .* c)) .^ (1 ./ 3);
templ = s + 1.0 ./ s + 1.0;
P = F ./ (3.0 .* templ .* templ .* G .* G);
Q = sqrt(1.0 + 2.0 .* e4 .* P);
r0 = -P .* e2 .* r ./ (1.0 + Q) + ...
    sqrt(abs(0.5 .* a2 .* (1.0 + 1.0 ./ Q) - ...
    P .* ome2 .* z2 ./ (Q .* (1.0 + Q)) - 0.5 .* P .* r2));
temp2 = r - e2 .* r0;
U = sqrt(temp2 .* temp2 + z2);
V = sqrt(temp2 .* temp2 + ome2 .* z2);
z0 = b2 .* z ./ (a .* V);

% calculate latitude
lat(valid) = atan2(z(valid) + e_b2 .* z0(valid), r(valid))*180/pi; % atan2d not available until MATLAB 2012b

% calculate altitude
alt(valid) = U(valid) .* (1.0 - b2 ./ (a .* V(valid)));

if nargout < 2 % Matrix/vector form, rather than componentwise form, was requested
    lat = [lat; lon; alt];
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////