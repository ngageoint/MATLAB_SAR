function ecf_value = ned_to_ecf(ned_value, orp_ecf, position_boolean)
%NED_TO_ECF Convert North-East-Down coordinates to ECF 
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Convert north, east, and down coordinates to ECF (Earth Centered Fixed).
%
% USAGE:
%   ecf_value = ned_ecf(ned_value, orp_ecf)
%
% INPUTS:
%   ned_value - required : north-east-down coordinates                [m, m, m]
%   orp_ecf - required : ecf coordinates of ORP                       [m, m, m]
%   position_boolean - optional : Is ecf_value a position coordinate? [true, false]
%      Unit normal vectors, velocity vectors, etc. do not require ORP offset.
%
% OUTPUTS:
%   ecf_value - required : ecf coordinates of ned_value               [m, m, m]

ecf_value=(ecf_ned_rot_mat(orp_ecf).')*ned_value(:);
if nargin<3||position_boolean % Assume value is a position coordinate, unless told otherwise
    ecf_value=ecf_value+orp_ecf(:);
end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////