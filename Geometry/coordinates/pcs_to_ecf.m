function ecf_value = pcs_to_ecf(pcs_value, orp_ecf, x_axis_orientation, position_boolean)
%PCS_TO_ECF Convert Processor Coordinates System to ECF 
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Convert PCS local coordinates to ECF (Earth Centered Fixed).
%
% USAGE:
%   ecf_value = pcs_to_ecf(ned_value, orp_ecf, x_axis_orientation, position_boolean)
%
% INPUTS:
%   pcs_value - required : PCS coordinates                            [m, m, m]
%   orp_ecf - required : ecf coordinates of reference point           [m, m, m]
%   x_axis_orientation : PCS X-axis orientation  [degrees clockwise from north]
%   position_boolean - optional : Is ecf_value a position coordinate? [true, false]
%      Unit normal vectors, velocity vectors, etc. do not require ORP offset.
%
% OUTPUTS:
%   ecf_value - required : ecf coordinates of ned_value               [m, m, m]

if isvector(pcs_value) % Assure correct orientation
    pcs_value = pcs_value(:);
end

% Calculate intermediates
orp_lla = ecf_to_geodetic(orp_ecf);
cos_lat = cosd(orp_lla(1));
sin_lat = sind(orp_lla(1));
cos_lon = cosd(orp_lla(2));
sin_lon = sind(orp_lla(2));
cos_phi = cosd(x_axis_orientation);
sin_phi = sind(x_axis_orientation);

% Transformation matrix
T = [-(sin_lat*cos_lon*cos_phi)-(sin_lon*sin_phi),...
        -(sin_lat*cos_lon*sin_phi)+(sin_lon*cos_phi),...
        cos_lat*cos_lon;...
    -(sin_lat*sin_lon*cos_phi)+(cos_lon*sin_phi),...
        -(sin_lat*sin_lon*sin_phi)-(cos_lon*cos_phi),...
        cos_lat*sin_lon;...
    cos_lat*cos_phi,...
        cos_lat*sin_phi,...
        sin_lat];
ecf_value = T*pcs_value;

if nargin<4||position_boolean % Assume value is a position coordinate, unless told otherwise
    ecf_value=ecf_value+repmat(orp_ecf(:),1,size(pcs_value,2));
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////