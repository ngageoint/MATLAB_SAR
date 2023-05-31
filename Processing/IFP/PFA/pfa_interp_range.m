function [phase_history] = pfa_interp_range(phase_history, k_a, k_r0, k_r_ss, new_v, interp_type)
%PFA_INTERP_RANGE Range interpolation stage of the polar formatting
%algorithm
%
% Inputs:
%     phase_history: Complex pulse sample data, motion-compensated phase
%                    history.  The phase history is expected to be supplied
%                    with each pulse as a column.
%     k_a:           k_a is the angle between each pulse and a reference
%                    position.  k_a could be viewed as the "squint" angle
%                    for each pulse for a sensor flying a straight-line.
%     k_r0:          Start radial position of each pulse.  Can be either a
%                    scalar that is constant for all pulses, or a vector
%                    containing a start frequency for each pulse.
%     k_r_ss:        Radial position spacings between the samples in each
%                    pulse.  Can be either a scalar that is constant for
%                    all pulses, or a vector containing a frequency step
%                    size for each pulse.
%     new_v:         Vector containing the V coordinates to interpolate to.
%
% We view the pulse sample points in a U,V coordinate system in which the V
% axis is parallel to the center pulse and directed 'outward' (increasing
% frequency) in the k_x, k_y plane.  U is orthoganal to V and increases to
% the left (when viewed from above).  We will make the units of our k-space
% to be cycles/meter, which is consistent with SICD metadata.
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

% Expand possibly scalar arguments
k_r0 = k_r0.*ones(size(phase_history,2),1);
k_r_ss = k_r_ss.*ones(size(phase_history,2),1);

for pulse=1:size(phase_history,2)
    % Compute the k_r values.  k_r for each sample is the radial position
    % or the distance in k-space from the origin.
    k_r = k_r0(pulse) + (k_r_ss(pulse)*(0:(size(phase_history,1)-1)));

    % Compute the k_v coordinates.  This is simply the rectangular
    % coordinate of the polar k_r, k_a coordinates in the V direction.
    k_v = k_r * cos(k_a(pulse));

    % We now have the phase history samples and the U,V coordinates of
    % each sample.  We now need to interpolate the data to a
    % rectangular grid.
    if strcmpi(interp_type,'sinc')
        phase_history(:,pulse) = sinc_interp(phase_history(:,pulse).', k_v, new_v.');
    else
        phase_history(:,pulse) = interp1(k_v, phase_history(:,pulse), new_v, interp_type, 0.0);
    end    
end

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////