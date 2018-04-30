function [ resolution, extent, delta_azimuth, total_azimuth ] = ...
    pulse_info_to_resolution_extent( range_vectors, center_frequency, delta_frequency, bandwidth, num_pulses )
%PULSE_INFO_TO_RES_EXTENT Computes theoretical resolution and maximum scene
%size for a set of SAR pulses
%
% [resolution, extent] = pulse_info_to_resolution_extent( range_vectors,...
% center_frequency, delta_frequency, bandwidth )
%
% Takes an array of RANGE_VECTORS for each pulse in meters (individual
% range vectors are stored x,y,z in rows).  Also requires CENTER_FREQUENCY
% (Hz), DELTA_FREQUENCY (step size of frequency data (Hz)), and BANDWIDTH
% (Hz), to compute theoretical scene RESOLUTION ([range cross-range],
% null-to-null in meters) and maximum scence EXTENT ([range cross-range],
% in meters). CENTER_FREQUENCY, DELTA_FREQUENCY, and BANDWIDTH can be a
% vector containing values for each pulse or single values that hold for
% all pulses.
%
% Values for resolution and extent are only approximate.  Does not project
% range vectors into an image plane or consider polar inscription.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

c = 299792458; % Speed of light (m/s)

% One could pass per-pulse metadata for every pulse, or just pass the first
% and last pulses, in which case we need to know the number of pulses.
if nargin<5, num_pulses = length(range_vectors); end

% Compute unit range vectors
range_vectors_norm=repmat(sqrt(sum(range_vectors.^2,2)),[1 3]);
unit_range_vectors=range_vectors./range_vectors_norm;

% Determine the total azimuth angle of the aperture (radians)
% Assumes a total azimuth less than 180 degrees
total_azimuth=acos(sum(unit_range_vectors(1,:).*unit_range_vectors(end,:),2));

% Determine the average azimuth angle step size (radians)
delta_azimuth=total_azimuth/num_pulses; % Approximation assuming pulses evenly spaced in angle
% Another way to do it if all pulse metadata is given:
% delta_azimuth=acos(mean(sum(unit_range_vectors(1:end-1,:).*unit_range_vectors(2:end,:),2)));

% Determine the maximum scene size of the image (m)
extent = [c/(2*mean(delta_frequency))... % Range
          c/(2*delta_azimuth*mean(center_frequency))]; % Cross-range

% Determine the resolution of the image (m)
resolution = [c/(2*mean(bandwidth))... % Range
              c/(2*total_azimuth*mean(center_frequency))]; % Cross-range

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////