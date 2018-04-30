function rot_mat = ecf_ned_rot_mat(orp_ecf)
%ECF_NED_ROT_MAT Calculate ECF to North-East-Down rotation matrix
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Calculate ECF (Earth Centered Fixed) to north-east-down rotation matrix.
%
% USAGE:
%   rot_mat = ecf_ned_rot_mat(orp_ecf)
%
% INPUTS:
%   orp_ecf - required : ecf coordinates of ORP                       [m, m, m]
%
% OUTPUTS:
%   rot_mat - required : rotation matrix

orp_lla=ecf_to_geodetic(orp_ecf);
rot_mat=[cosd(-90-orp_lla(1)), 0, -sind(-90-orp_lla(1));...
    0, 1, 0;...
    sind(-90-orp_lla(1)), 0, cosd(-90-orp_lla(1))]*...
    [cosd(orp_lla(2)), sind(orp_lla(2)), 0;...
    -sind(orp_lla(2)), cosd(orp_lla(2)), 0;...
    0, 0, 1];

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////