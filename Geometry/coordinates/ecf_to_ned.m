function ned_value = ecf_to_ned(ecf_value, orp_ecf, position_boolean)
%ECF_TO_NED Convert ECF coordinates to North-East-Down
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////
%
% Convert ECF (Earth Centered Fixed) coordinates to north, east, and down.
%
% USAGE:
%   ned_value = ecf_to_ned(ecf_value, orp_ecf)
%
% INPUTS:
%   ecf_value - required : ecf x, y, z coordinates                    [m, m, m]
%   orp_ecf - required : ecf coordinates of ORP                       [m, m, m]
%   position_boolean - optional : Is ecf_value a position coordinate? [true, false]
%      Unit normal vectors, velocity vectors, etc. do not require ORP offset.
%
% OUTPUTS:
%   ned_value - required : north-east-down coordinates of ecf_value   [m, m, m]

if nargin<3||position_boolean % Assume value is a position coordinate, unless told otherwise
    ecf_value=ecf_value-orp_ecf(:);
end
ned_value=ecf_ned_rot_mat(orp_ecf)*ecf_value;

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////