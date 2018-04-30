function [ degrees, minutes, seconds] = latlonvec(decimal_degrees)
%LATLONVEC Convert decimal degrees into degrees/minutes/seconds
%
% Description: 
% Return degrees, minutes, seconds given decimal degrees. Sign of
% coordinate will be returned in degrees portion of coordinate.
%
% Sample fuction call: 
% [ degrees, minutes, seconds] = latlonvec(decimal_degrees)
%
% Author: Scott Lee
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

degrees = fix(abs(decimal_degrees));	% integer degrees
frac = abs(decimal_degrees) - degrees;	% fractional degrees, used to compute minutes 
minutes = fix(frac*60);	    % integer minutes
frac = frac - minutes/60;	% fractional minutes, used to compute seconds
seconds = frac*3600;		% decimal seconds

% Handle sign.  Degrees portion will contain the sign of the coordinate.
% Minutes and seconds will always be positive.
% sign function returns -1, 0, +1 for x < 0, x == 0, x > 0, respectively
degrees = sign(decimal_degrees).*degrees;

if nargout<2
  degrees = [degrees minutes seconds];
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////