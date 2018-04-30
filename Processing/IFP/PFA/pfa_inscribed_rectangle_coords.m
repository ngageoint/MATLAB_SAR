function [k_v_bounds, k_u_bounds] = pfa_inscribed_rectangle_coords(k_a, k_r0, bw)
%PFA_INSCRIBED_RECTANGLE_COORDS Computes the coordinates for a rectangle
% that inscribes the SAR polar annulus.
%
% This function assumes that the data has both postive and negative polar
% angles, and that the inscribed rectangle is oriented with the U/V
% directions.  Also assumes all functions have roughly similar frequency
% ranges (won't work for "step chirp" data).
%
% Inputs:
%     k_a:            k_a is the polar angle between each pulse and a
%                     reference position.
%     k_r0:           Start frequency (or minimum radial position) of each
%                     pulse.  Can be either a scalar that is constant for
%                     all pulses, or a vector containing a start frequency
%                     for each pulse.
%     bw:             Total effective bandwidth (or radial extent in
%                     "K-space" in the image plane) of each pulse.  Can be
%                     either a scalar that is constant for all pulses, or a
%                     vector containing bandwidths for each pulse.
%
% Outputs:
%     k_v_bounds:     Upper and lower limit of the inscribed rectangle in
%                     the V dimension (top and bottom edges).
%     k_u_bounds:     Upper and lower limit of the inscribed rectangle in
%                     the U dimension (left and right edges).
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
    
%% Compute the bounds of the inscribed rectangle in k-space.
% The bottom of the inscribed rectangle.
k_v_bounds = max(k_r0 .* cos(k_a));
% The sides, where the bottom of the inscribed rectangle bumps into angular
% bounds defined by the first/last pulses.
k_u_bounds = k_v_bounds * tan([min(k_a) max(k_a)]);
% "Radial" extent of the samples in each pulse in cycles per meter.
k_r_high = k_r0 + bw;
% If start frequencies and bandwidths were exactly the same for every
% pulse, the top of the inscribed rectangle could be solved analytically:
% k_v_bounds(2) = min(sqrt((min(k_r_high)^2)-(k_u_bounds.^2)));
% "Keystone" datasets can be simply computed this way:
% k_v_bounds(2) = min(k_r_high .* cos(k_a));
% However, a more generic solution that includes both of the above cases
% and datasets whose start frequencies and/or bandwidths can vary slightly
% per pulse is used here:
% u/v positions of the highest radial frequency each pulse
k_u_high = k_r_high .* sin(k_a);
k_v_high = k_r_high .* cos(k_a);
% Pulses whose highest radial frequency falls within k_u_bounds
valid = (k_u_high>=k_u_bounds(1)) & (k_u_high<=k_u_bounds(2));
% Top of the box can't be above any of these v frequencies
k_v_bounds(2) = min(k_v_high(valid));

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////