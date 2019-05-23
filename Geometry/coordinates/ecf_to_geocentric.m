function pos_lla = ecf_to_geocentric(pos_ecf)
%ECF_TO_GEOCENTRIC Convert ECF coordinates to geocentric
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Convert ECF (Earth Centered Fixed) coordinates to geocentric latitude, 
% longitude, and altitude.
%
% USAGE:
%   pos_lla = ecf_to_geocentric(pos_ecf)
%
% INPUTS:
%   pos_ecf - required : ecf x, y, z coordinates [m, m, m]
%
% OUTPUTS:
%   pos_lla - required : geocentric latitude, longitude, and altitude [deg, deg, m]
%
% NOTES:
%   Bate, Roger. Mueller, Donald. White, Jerry. Fundamentals of Astrodynamics. 
%   Dover Publications, Inc. 1971.
%
% VERSION:
%   1.0
%     - Sean Hatch 20070911
%     - initial version
%
% TODO:

% define constants
e2 = 6.6943799901377997e-3;  % eccentricity squared of Earth (WGS 84 value)
a = 6378137.0;               % semimajor radius of the Earth (WGS 84 value)
rad_to_deg = 180 ./ pi;

% calculate derived constants
ome2 = 1.0 - e2;
a2 = a .* a;
b = a .* sqrt(ome2);         % semiminor radius of the Earth
b2 = b .* b;

% calculate intermediates
x = pos_ecf(1);
y = pos_ecf(2);
z = pos_ecf(3);
x2 = x .* x;
y2 = y .* y;
z2 = z .* z;
r2 = x2 + y2;
r = sqrt(r2);

% handle valid solution
if ((r2 > 0.0) || (z2 > 0.0))
    
    % calculate longitude
    lon = atan2(y, x);

    % calculate latitude
    lat = atan2(z, r);

    % calculate distance to surface of ellipsoid
    R = sqrt(1.0 ./ (r2 ./ a2 + z2 ./ b2)) .* sqrt(r2 + z2);

    % calculate altitude
    alt = sqrt(r2 + z2) - R;

% handle invalid solution
else
    
    % assign default values
    lon = 0.0;
    lat = 0.0;
    alt = 0.0;
    
end

% convert latitude and longitude to degrees
lat = lat .* rad_to_deg;
lon = lon .* rad_to_deg;

% prepare the result
pos_lla = [lat; lon; alt];

% all done
return;

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////