function [phase_history] = pfa_interp_azimuth(phase_history, k_a, v_coords, new_u, interp_type)
%PFA_INTERP_AZIMUTH Azimuth interpolation stage of the polar formatting algorithm
%
% Inputs:
%     phase_history: Trapezoidal phase history data that has common
%                    coordinates in the V dimension across all pulses
%                    (either because it was interpolated that way, or
%                    because it was collected that way.)  Each pulse is a
%                    column.
%     k_a:           k_a is the angle between each pulse and a reference
%                    pulse.  k_a could be viewed as the "squint" angle for
%                    each pulse for a sensor flying a straight-line.
%     v_coords:      Vector containing the V coordinates for each sample.
%                    Since are input phase history is trapezoidal, this
%                    should be the same for every pulse.
%     new_u:         Vector containing the U coordinates to interpolate to.
%
% We view the pulse sample points in a U,V coordinate system in which the V
% axis is parallel to the center pulse and directed 'outward' (increasing
% frequency) in the k_x, k_y plane.  U is orthoganal to V and points in the
% direction of increasing pulse number-- like that on page 185 of Carrara,
% Goodman, Majewski.  We will make the units of our k-space to be
% cycles/meter, which is consistent with SICD metadata.
%
% In this function, (k_r, k_a) are the polar coordinates (radial and
% angular) of every pulse/sample, and (k_u, k_v) are the rectangular
% coordinates in the U,V coordinate system.
%
% Authors: Wade Schwartzkopf and Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

if ~exist('interp_type','var')
    interp_type = 'linear'; % Ideally would be sinc
end

for sample = 1:size(phase_history,2)
    % Compute the k_u rectangular coordinates.  When we interpolated
    % radially to a known V coordinate (v_coords), that also resulted in
    % updated U coordinates.  We compute those new U coordinates here.
    k_u = tan(k_a) .* v_coords(sample);

    % For each sample location interpolate to new desired values of U.
    if strcmpi(interp_type,'sinc')
        phase_history(:,sample) = sinc_interp(phase_history(:,sample).', k_u.', new_u.');
    else
        phase_history(:,sample) = interp1(k_u, phase_history(:,sample), new_u, interp_type, 0.0);
    end
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////