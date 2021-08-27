function [data_block] = pfa_mem(phase_history, nbdata, sample_rate)
%PFA_MEM Implements the polar format IFP algorithm on an array of phase
%history data
%
% Inputs:
%     phase_history: Complex pulse sample data as returned as the first
%                    parameter of the read_cphd method of an open_ph_reader
%                    object.  The phase history is expected to be supplied
%                    with each pulse as a column.
%     nbdata:        Per-pulse narrowband data.  This is the same structure
%                    returned as the second return parameter of the
%                    read_cphd method of an open_ph_reader object.
%     sample_rate:   Samples per IPR.  Default is 1.5.
%
% We'll generate the pulse sample points in a U,V coordinate system in
% which the V axis is parallel to the center of aperture and directed
% 'outward' (increasing frequency) in the k_x, k_y plane.  U is orthoganal
% to V and increases to the left (when viewed from above).
%
% In this function, (k_r, k_a) are the polar coordinates (radial and
% angular) of every pulse/sample, and (k_u, k_v) are the rectangular
% coordinates in the U,V coordinate system.  We will make the units of our
% k-space to be cycles/meter, which is consistent with SICD metadata.
%
% Authors: Wade Schwartzkopf and Tom Krauss, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Parse input parameters
if nargin<3
    sample_rate = 1.5;
end

if any(any(diff(nbdata.SRPPos)))  % Assure spotlight data
    error('PFA_MEM:UNSUPPORTED_COLLECT_TYPE','Unsupported collection mode.  Currently only spotlight data is supported.');
else
    scp = nbdata.SRPPos(1,:);
end

[ bi_pos, bi_freq_scale ] = pfa_bistatic_pos( nbdata.TxPos, nbdata.RcvPos, nbdata.SRPPos );
% Temporary for testing fpn/ipn.  Actually this should be an input argument.
fpn=wgs_84_norm(scp).'; % Compute normal to WGS_84 ellipsoid
scp = [0 0 0]; % bistatic position is with respect to scp as origin
ref_pulse = ceil(size(phase_history,2)/2); % Center pulse seems as good as any...
arp_coa_vel = diff(bi_pos(ref_pulse + [1 -1],:));
arp_coa = bi_pos(ref_pulse,:);
srv=(arp_coa-scp).'; % slant range vector
look=fpn*cross(srv,arp_coa_vel).';
ipn=look*cross(srv,arp_coa_vel); ipn=ipn/norm(ipn); % Slant plane unit normal
% end temporary portion here

% Compute angular position of each pulse
[k_a, k_sf] = pfa_polar_coords(bi_pos, scp, arp_coa, ipn, fpn); % Angular coordinate of each pulse
% Converting from raw RF frequency of the received pulse to radial position
% of each pulse/sample in the image formation plane in "K-space" involves a
% few scaling factors:
rf_to_rad = (2/SPEED_OF_LIGHT) .* ... % Convert from cycles/second to cycles/meter
    k_sf .* ... % Compensate for out-of-plane motion by projecting into image formation plane
    bi_freq_scale; % For bistatic collects, a factor to account for using the equivalent monostatic position
k_r0 = nbdata.SC0 .* rf_to_rad; % Radial position of the first sample in each pulse
k_r_ss = nbdata.SCSS .* rf_to_rad; % Radial position spacings between the samples in each pulse
% Compute new coordinates onto which to interpolate
[k_v_bounds, k_u_bounds] = pfa_inscribed_rectangle_coords(k_a, k_r0, ...
    k_r_ss * (size(phase_history,1)-1));
% Coordinates onto which we will interpolate
new_v = linspace(k_v_bounds(1), k_v_bounds(2), size(phase_history,1)).';
new_u = linspace(k_u_bounds(1), k_u_bounds(2), size(phase_history,2)).';

% The interpolation will be done in two parts.  First we'll interpolate
% radially within each pulse to a constant spacing in V, then we'll use the
% just-interpolated values to interpolate across pulses (at each sample) to
% a constant spacing in U.
data_block = pfa_interp_range(phase_history, k_a, k_r0, k_r_ss, new_v);
data_block = pfa_interp_azimuth(data_block.', k_a, new_v, new_u);

% FFT 
data_block = pfa_fft_zeropad_1d(data_block.', sample_rate); % FFT V
data_block = pfa_fft_zeropad_1d(data_block.', sample_rate); % FFT U

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////