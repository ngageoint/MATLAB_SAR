function gotcha_pr2cphd( directory_name, output_file, pass, pol, azimuth_range )
%GOTCHA_PR2CPHD Converts a dataset from the AFRL 2D/3D SAR Volumetric
%public release dataset to CPHD file format.
%
% Usage: gotcha_pr2cphd( directory_name, output_file, pass, pol, azimuth_range)
%
% Inputs: directory_name: Data root directory for the GOTCHA public release
%                         dataset.  Assumes all files are kept in the same
%                         directory structure that the dataset is
%                         distributed with.
%         output_file:    Name of CPHD file to save.
%         pass:           What pass to image (1-8).  Default = 1.
%         pol:            What polarization to image ('HH','HV','VH','VV').
%                         Default = 'HH'.
%         azimuth_range:  Two-element vector of the first and last azimuth
%                         angle (in degrees) to convert to CPHD format.
%                         Default = [1 360] (all available angles).
%
% This function is just a wrapper to CONVERT_TO_CPHD, but allows for
% simpler selection of the pass, polarization, and azimuth range desired.
%
% Author: Wade Schwartzkopf, NGA/IDT
%
% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////

% Setup default parameters
if nargin<3
    pass = 1;
end
if nargin<4
    pol = 'HH';
end
if nargin<5
    azimuth_range = [1 360];
end

% Calculate which pulse numbers go with which files in directory for
% desired azimuth range.
last_pulse = 0;
pulse_ranges = zeros(azimuth_range(2),2); % First and last pulse number in every file
for az = 1:azimuth_range(end)
    fname = fullfile(directory_name, sprintf('pass%d',pass),pol,...
        sprintf('data_3dsar_pass%d_az%03d_%s.mat',pass,az,pol));
    load(fname); % Loads data into a structure named 'data'
    num_pulses = size(data.fp,2);
    pulse_ranges(az,:) = last_pulse + [1 num_pulses];
    last_pulse = last_pulse + num_pulses;
end
first_pulse = pulse_ranges(azimuth_range(1),1);
last_pulse = pulse_ranges(azimuth_range(2),2);

convert_to_cphdx(fname,output_file,'pulse_indices',first_pulse:last_pulse);
% convert_to_cphd30(fname,output_file,first_pulse:last_pulse); % Old way

end

% //////////////////////////////////////////
% /// CLASSIFICATION: UNCLASSIFIED       ///
% //////////////////////////////////////////